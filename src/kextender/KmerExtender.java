/*
 * Copyright 2015 Australian Centre For Plant Functional Genomics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package kextender;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.io.BufferedWriter;
import shared.Reporter;
import shared.InputReaderProducer;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentNavigableMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import shared.StdRedirect;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class KmerExtender {

//    private PairMersMap pairMersMap = new PairMersMap(); //TODO user to give expected num of elems?
    private final ArrayList<String> inputFileNamesList = new ArrayList<>();
//    private Integer KMER_LENGTH;
    private Integer KMER_LENGTH_MIN;
    private Integer KMER_LENGTH_MAX;
    private Integer KMER_LENGTH_STEP;
    private SeedSequences seedSequences;
    private Integer MIN_KMER_FREQUENCY = 1; //PLANNED but not really used other than to exclude , must be >=1
    private int MAX_THREADS;
    private int INPUT_BUFFER_SIZE;
    private int INPUT_QUEUE_SIZE;

    private boolean OUTPUT_FASTA = false;
    private String NAME_PREFIX = "";
//    private String DEBUG_FILE;
//    private String STATS_FILE;
//    private Integer HASH_ARRAY_SIZE = Integer.MAX_VALUE - 3;
//    private Integer MULTIPASS_COMPRESS = null;

    private final String TOOL_NAME;
    private final int HELP_WIDTH = 140;
    private int WRITER_BUFFER_SIZE = 8192;

    private enum InputType {

        KMERS, FASTA, FASTQ;
    }

    public KmerExtender(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        readArgValues(optSet);
//        if (RUN_SOME_WILD_AND_WONDERFUL_STUFF) {
//            new kmerextender.ideas.Alternative(MAX_THREADS, inputFileNamesList, KMER_LENGTH);
//        } else {
        runKmerExtender(optSet);
//        }
    }

    private void readArgValues(OptSet optSet) {
        new StdRedirect(optSet, TOOL_NAME, StdRedirect.RedirectType.REDIRECT_ERR);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        INPUT_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        INPUT_QUEUE_SIZE = (int) optSet.getOpt("Q").getValueOrDefault();

//        MIN_KMER_FREQUENCY = (int) optSet.getOpt("min-frequency").getValueOrDefault();
        if (optSet.getOpt("k").isUsed()) {
            setKmerLength((int) optSet.getOpt("k").getValueOrDefault());
        } else {
            if (optSet.getOpt("k-mer-min").isUsed()) {
                KMER_LENGTH_MIN = (int) optSet.getOpt("k-mer-min").getValueOrDefault();
            }
            if (optSet.getOpt("k-mer-max").isUsed()) {
                KMER_LENGTH_MAX = (int) optSet.getOpt("k-mer-max").getValueOrDefault();
            }
            KMER_LENGTH_STEP = (int) optSet.getOpt("k-mer-step").getValueOrDefault();
            if (KMER_LENGTH_MIN != null & KMER_LENGTH_MAX != null) {
                if (KMER_LENGTH_MAX < KMER_LENGTH_MIN) {
                    Reporter.report("[ERROR]", "The k-mer-min value (" + KMER_LENGTH_MIN + ") set higher than k-mer-max value (" + KMER_LENGTH_MAX + "), adjusting to " + KMER_LENGTH_MAX + ",", TOOL_NAME);
                    KMER_LENGTH_MAX = KMER_LENGTH_MIN;
                }
            } else if (KMER_LENGTH_MIN == null) {
                KMER_LENGTH_MIN = KMER_LENGTH_MAX;
            } else if (KMER_LENGTH_MAX == null) {
                KMER_LENGTH_MAX = KMER_LENGTH_MIN;
            }
        }

        if (optSet.getOpt("seed-file").isUsed()) {
            seedSequences = new SeedSequences((String) optSet.getOpt("seed-file").getValueOrDefault(), optSet.getOpt("purged-input").isUsed());
            if (seedSequences.getSeedSequences().isEmpty()) {
                Reporter.report("[ERROR]", "No seeds found in " + (String) optSet.getOpt("seed-file").getValueOrDefault() + ", proceeding without...", TOOL_NAME);
                if (KMER_LENGTH_MIN != KMER_LENGTH_MAX) {
                    Reporter.report("[FATAL]", "Varied k implemented  (for now) for seed extension only - no seed(s) found", TOOL_NAME);
                    System.exit(1);
                }
            }
        }

//        MIN_KMER_FREQUENCY = (int) optSet.getOpt("k").getValueOrDefault();
        MAX_THREADS = (int) optSet.getOpt("t").getValueOrDefault();

        if (optSet.getOpt("f").isUsed()) {
            OUTPUT_FASTA = true;
        }
        if (optSet.getOpt("p").getValueOrDefault() != null) {
            NAME_PREFIX = (String) optSet.getOpt("p").getValueOrDefault();
        }

//        String outRedirect;
//        String errRedirect;
//        if ((outRedirect = (String) optSet.getOpt("o").getValueOrDefault()) != null) {
//            try {
//                File file = new File(outRedirect);
//                PrintStream printStream;
//                printStream = new PrintStream(new FileOutputStream(file));
//                System.setOut(printStream);
//            } catch (FileNotFoundException ex) {
//                Reporter.report("[ERROR]", "Failed redirecting stdout to " + outRedirect, TOOL_NAME);
//            }
//        }
//        if ((errRedirect = (String) optSet.getOpt("e").getValueOrDefault()) != null) {
//            try {
//                File file = new File(errRedirect);
//                PrintStream printStream;
//                printStream = new PrintStream(new FileOutputStream(file));
//                System.setErr(printStream);
//            } catch (FileNotFoundException ex) {
//                Reporter.report("[ERROR]", "Failed redirecting stderr to " + errRedirect, TOOL_NAME);
//            }
//        }
//        for(Opt o: optSet.getOptsList()) {
//            Reporter.report("[INFO]", o.getOptLabelString()+" "+o.getValueOrDefault(), toolName);
//        }
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFileNamesList.addAll(po.getValues());
            }
        }
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("(k)extender is designed for rapid identification of unambiguous extensions (a.k.a. unitigs) "
                + "from a set of k-mers. Note that due to specifics of the implementation there is no advantage in using odd k values. "
                + "Frequncies of k-mers are not currently being utilized, so it is best to use input which is unlikely to contain erroneous k-mers. "
                + "");
//                + "Each input k-mer is split and stored twice (as k-1-mer and the remaining on-base extension). "

        //INPUT
        optSet.setListingGroupLabel("[Input settings - general extender]");
        optSet.addOpt(new Opt('k', "k-mer-length", "Required only if input other than a list of k-mers", 1).setMinValue(4).setMaxValue(2048));
//        optSet.addOpt(new Opt('M', "min-frequency", "[IMPLEMENTING]", 1).setMinValue(1).setMaxValue((int)Byte.MAX_VALUE).setDefaultValue(1));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of records (k-mers or FASTQ reads or pairs depending on input) "
                + "passed to in-queue", 1024, 1, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for populator threads to pick-up",
                2, 1, 256));
//        optSet.addOpt(new Opt('m', "min-frequency", "....", 1,1,Integer.MAX_VALUE));
        int footId = 1;
        String foot = "Warning! Exploring a large range of k values for a significant input [k-mers/FAST[A|Q]] "
                + "will make the memory requirements explode, "
                + "as it is done in parallel to avoid excessive I/O and to preserve stdin handling";
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Variable k-mer size for longest extension of a seed]");
        optSet.addOpt(new Opt('S', "seed-file", "Fasta file containing \"seeds\" to be extended", 1));
//        optSet.addOpt(new Opt(null, "exclude", "Blacklisted k-mers to be excluded before extending", 1));
        optSet.addOpt(new Opt(null, "purged-input", "The input k-mers do not contain any of the seed-mers"));
        optSet.addOpt(new Opt(null, "k-mer-min", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-max", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-step", "", 1).setMinValue(1).setDefaultValue(2).addFootnote(footId, foot));
//        optSet.addOpt(new Opt(null, "k-values", "Explicitly set k values to be explored", 2).setMinValue(1).setDefaultValue(2).addFootnote(footId, foot));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        int procs = Runtime.getRuntime().availableProcessors();
        optSet.addOpt(new Opt('t', "threads", "Number of worker threads, at various stages of program execution additional 1-3 threads may be running",
                procs, 1, Math.min(procs, 255))); //max 255 threads as w use a byte to store thread id on visited pairmers with Byte.MAX_VALUE reserved as intial value representing not visited

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('m', "min-length", "Do not output extended sequnces shorter than <arg>", 1).setMinValue(5).setDefaultValueDescription("k+1"));
        optSet.addOpt(new Opt('f', "fasta-out", "Output each k-mer as a separate FASTA record instead of just listing extended nucleotide sequences"));
        optSet.addOpt(new Opt('p', "fasta-id-prefix", "Prefix each FASTA identifier with <arg> ", 1));
        optSet.addOpt(new Opt('o', "out-file", "Print extended sequences to <arg> file", 1).setDefaultValue("/dev/stdout"));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
//        optSet.addOpt(new Opt(null, "ambiguous-only", "Only output extensions around ambiguous positions", 1));
//        optSet.addOpt(new Opt(null, "blunt-only", "Only output  unambiguous extensions (ones that end due to lack of sequence information rather than alternatives)", 1));

        optSet.addOpt(new Opt('s', "stats-file", "Write extension stats to this file [TO BE RE-IMPLEMENTED - not much output here]", 1));
        optSet.addOpt(new Opt('d', "debug-file", "Write unkosher extensions details to this file - in practice these are just lists of k-mers for the unextended palindromic/circular sequences", 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    /**
     * Wrapper method executing individual steps
     */
    private void runKmerExtender(OptSet optSet) {
        Reporter.report("[INFO]", "Initialized, will use " + MAX_THREADS + " thread(s) to populate map ", TOOL_NAME);

        ArrayList<Integer> kSizes = new ArrayList<>();
        if (KMER_LENGTH_MIN != null) {
            for (int k = KMER_LENGTH_MIN; k < KMER_LENGTH_MAX + 1; k += KMER_LENGTH_STEP) {
                kSizes.add(k);
            }
        }
        Collections.sort(kSizes, Collections.reverseOrder());

        //READ k-mers AND POPULATE A MAP FOR EACH SIZE OF k
        PairMerMaps pairMerMaps = new PairMerMaps(kSizes, TOOL_NAME);
        readKmersAndPopulatePairMersMaps(pairMerMaps);

        Iterator<Integer> kIterator = kSizes.iterator();
        while (kIterator.hasNext()) {
            Integer k = kIterator.next();
//
//        }
//        for (Integer k : kSizes) {
            PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
//            Reporter.report("[INFO]", "Finished populating map, counting elements... ", TOOL_NAME);
            Reporter.report("[INFO]", "Finished populating map, k=" + k + ", n=" + NumberFormat.getNumberInstance().format(pairMersMap.size()), TOOL_NAME);
            if (pairMersMap.isEmpty()) {
                pairMerMaps.removePairMersMap(k);
                kIterator.remove();
            }
        }

        if (pairMerMaps.size() == 0) {
            Reporter.report("[WARN]", "Empty k-mer map(s), nothing to extend... terminating", TOOL_NAME);
            System.exit(0); //not sure if this should be non-zero 
        }

        //PURGING ALSO IDENTIFIES MOST TERMINAl PAIRMERS -> ALLOWING EFFICIENT MULTI-THREADED TRAVERSAL
        Reporter.report("[INFO]", "Purging map and identifying terminal PairMers...", TOOL_NAME);
        purgePopulatedPairMersMaps(pairMerMaps);

        ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers = null;
        if (seedSequences != null) {
            if (!optSet.getOpt("purged-input").isUsed()) {
                //PURGE NON-TERMINAL SeedPairMers FROM PairMerMaps
                purgeSeedMersFromPairMersMaps(kSizes, pairMerMaps);
            }
            //CREATE PairMers REPRESENTING SEED-ENDS
            Reporter.report("[INFO]", "Now populate seedMersMap", TOOL_NAME);
            kToSeedMers = populateSeedMersMaps(kSizes);
        }

        String outFile = (String) optSet.getOpt("out-file").getValueOrDefault();
        Reporter.report("[INFO]", "Now extending...", TOOL_NAME);
        //EXTEND - CAN BE PARALLELIZED IF NO INTENTION TO GENERATE CROSS-k-EXTENSIONS
        if (seedSequences == null) {
            for (Integer k : kSizes) {
                PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
                PairMersExtender pairMersExtender = new PairMersExtender((String) optSet.getOpt("debug-file").getValueOrDefault(),
                        (String) optSet.getOpt("stats-file").getValueOrDefault(), TOOL_NAME);
                Integer minLen = (int) optSet.getOpt("min-length").getValueOrDefault(k + 1);
                pairMersExtender.matchAndExtendKmers(k, pairMersMap, OUTPUT_FASTA, NAME_PREFIX,
                        outFile, minLen, MAX_THREADS);
                Reporter.report("[INFO]", "Finished extending for k=" + k, TOOL_NAME);

            }
//            pairMersExtender.matchAndExtendSeeds(k, pairMersMap, kToSeedMers.get(k));
        } else {
            extendSeedsParallel(kSizes, pairMerMaps, kToSeedMers, optSet);
        }

        //SELECT LONGEST EXTENSION
        if (seedSequences != null) {
            try {
                BufferedWriter out;
                if (outFile.endsWith(".gz")) {
                    out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile)), "UTF-8"), WRITER_BUFFER_SIZE);
                } else {
                    out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"), WRITER_BUFFER_SIZE);
                }
                for (SeedSequence seed : seedSequences.getSeedSequences()) {
//            Reporter.report("[INFO]", "Longest extension of the provided seed is from " + seed.getSequenceString().length() + " to " + longestSeedAfterExtension + " at k = " + kBest, TOOL_NAME);
//                Map.Entry<Integer, String> longest = seed.getLongestExtended();
                    Map.Entry<Integer, SeedExtensionsPair> longestExtensionLeft = seed.getLongestExtensionLeft(TOOL_NAME);
                    Map.Entry<Integer, SeedExtensionsPair> longestExtensionRight = seed.getLongestExtensionRight(TOOL_NAME);
//                StringBuilder output = new StringBuilder();
                    StringBuilder outputLR = new StringBuilder();
                    String extensionLeft = longestExtensionLeft.getValue().getExtensionLeft();
                    String extensionRight = longestExtensionRight.getValue().getExtensionRight();
                    if (OUTPUT_FASTA) {
                        outputLR.append(">");
                        if (!NAME_PREFIX.isEmpty()) {
                            outputLR.append(NAME_PREFIX).append("_");
                        }
                        outputLR.append(seed.getId());
                        outputLR.append(" L=").append(extensionLeft.length()).append(" at k=").append(longestExtensionLeft.getKey());
                        outputLR.append(" R=").append(extensionRight.length()).append(" at k=").append(longestExtensionRight.getKey());
                        outputLR.append(" extended=");
                        outputLR.append(extensionLeft.length() + seed.getSequenceString().length() + extensionRight.length());
                        outputLR.append(System.lineSeparator());
                        out.write(outputLR.toString());
                    }
//                output.append(longest.getValue());
//                System.out.println(output);
                    out.write(extensionLeft);
                    out.write(seed.getSequenceString());
                    out.write(extensionRight);
                    out.newLine();
//                    outputLR.append(extensionLeft).append(seed.getSequenceString()).append(extensionRight);
//                    System.out.println(outputLR);

                }
                out.close();
            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }
        }
        Reporter.report("[INFO]", "Finished extending k-mers", TOOL_NAME);
    }

    /**
     * Producer-consumer multi-threading to: - read the input (k-mers, FASTA,
     * etc) - generate PairMer objects - populate Map(s)
     */
    private void readKmersAndPopulatePairMersMaps(PairMerMaps pairMerMaps) {
        BlockingQueue inputQueue = new ArrayBlockingQueue(INPUT_QUEUE_SIZE);

        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = MAX_THREADS;
            Reporter.report("[INFO]", "Allocated " + threads + " thread(s) to map populating", TOOL_NAME);
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateExecutor = new ThreadPoolExecutor(threads + 1, threads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            //SPAWN INPUT READING THREAD
            ArrayList<Integer> kSizes = pairMerMaps.getkSizes();
            InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, kSizes, inputFileNamesList, INPUT_BUFFER_SIZE, TOOL_NAME);

            Future<?> future = readAndPopulateExecutor.submit(inputReaderProducer);
            futures.add(future);
            boolean splitInputSequenceIntoKmers = true;

            //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
            long timeStart = System.currentTimeMillis();
            int count = 0;
            while (inputReaderProducer.getGuessedInputFormat() == null) {
                try {
                    //IF nothing happens after 5 seconds
                    if (System.currentTimeMillis() - timeStart > 5000 && (count++ % 50 == 0)) {
                        Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", TOOL_NAME);
                    }
                    Thread.sleep(100); //wait for 1/10 of a second
                } catch (InterruptedException ex) {
                }
            }
//            Integer kSizeFromInput = null;
            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.InFormat.KMERS)) {
                splitInputSequenceIntoKmers = false;
            }

            //SPAWN THREADS TO POPULATE MAP
            for (int i = 0; i < threads; i++) {
                PairMerMapPopulatorConsumer consumer = new PairMerMapPopulatorConsumer(inputQueue, pairMerMaps,
                        splitInputSequenceIntoKmers, kSizes, MIN_KMER_FREQUENCY, false);
                futures.add(readAndPopulateExecutor.submit(consumer));
            }
            readAndPopulateExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                readAndPopulateExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!", TOOL_NAME);
            } catch (ExecutionException ex) {
                Reporter.report("[ERROR]", "PairMerSet populator execution exception!", TOOL_NAME);
                ex.printStackTrace();

            } catch (TimeoutException ex) {
                Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
                Reporter.report("[ERROR]", "PairMerSet populator timeout exception!", TOOL_NAME);
            }

//            if (pairMersMap.isOutOfMemory()) {
//                Reporter.report("[ERROR]", "Terminating, out of memory while populating the Map, k=" + KMER_LENGTH + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);
//                System.exit(1);
//            }
        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error!", TOOL_NAME);
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            System.exit(1);
        }

//        return kReturn;
    }

    private void purgePopulatedPairMersMaps(PairMerMaps pairMerMaps) {
        //PREPARE EXECUTOR SERVICE
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService purgeMapsExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

//        //No NEED TO SPLIT MAP(s) FOR PURGING IF #MAPS >= THREADS
//        int subMapThreads = pairMerMaps.size() >= threads ? 1 : threads;
        //INIT PURGER THREADS
        ArrayBlockingQueue<PairMersMap> mapsQueue = new ArrayBlockingQueue(Math.max(pairMerMaps.size(), threads) + 1);
        for (int i = 0; i < threads; i++) {
            futures.add(purgeMapsExecutorService.submit(new PairMerMapPurgerConsumer(mapsQueue, TOOL_NAME, MIN_KMER_FREQUENCY)));
        }

        //STATS OF MAP SPLITTING FOR MULTITHREADED PURGING
//        long total = 0;
//        long min = Long.MAX_VALUE;
//        long max = Long.MIN_VALUE;
//        int[] sizes = null;
//        ArrayList<ConcurrentNavigableMap<PairMer, PairMer>> mapChunks = new ArrayList<>();
        int BUFFER_SIZE = 100000;

        //PREPARE MAPS ON QUEUE TO BE PICKED-UP BY THREADS
        Iterator<Integer> it = pairMerMaps.getkSizes().iterator();
        try {
            while (it.hasNext()) {
                Integer k = it.next();
                PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
                if (pairMersMap.storeSize() == 0) { //Record for stats after purging, also continue if empty
                    continue;
                }
                if (pairMerMaps.size() >= threads) {
                    mapsQueue.put(pairMersMap);
                } else { //WE HAVE MORE THREADS THAN MAPS SO LETS SPLIT THE MAPS
//                    Reporter.report("[EXPERIMENTAL]", "Running recursive map splitting, need boundary cases testing as it can get stuck or crash", TOOL_NAME);
//                    long attempts = pairMersMap.recursiveSplitMap(mapChunks, pairMersMap.size() / threads / 100, pairMersMap.size() / threads / 2, mapsQueue);
//                    Reporter.report("[EXPERIMENTAL]", "Finished running recursive map splitting, total attempts = " + attempts, TOOL_NAME);
                    Iterator<PairMer> pairMersIterator = pairMersMap.iterator();
                    PairMer head = null;
                    PairMer tail = null;
                    long counter = 0;
                    if (pairMersIterator.hasNext()) {
                        head = pairMersIterator.next();
                        counter++;
                    }
                    while (pairMersIterator.hasNext()) {
                        if (counter++ % BUFFER_SIZE == 0) {
                            tail = pairMersIterator.next();
                            ConcurrentNavigableMap subMap = pairMersMap.getPairMersMap().subMap(head, tail);
                            mapsQueue.put(new PairMersMap(k, subMap, pairMersMap));
                            head = tail;
                        } else {
                            pairMersIterator.next();
                        }
                    }
                    ConcurrentNavigableMap tailMap = pairMersMap.getPairMersMap().tailMap(head);
                    mapsQueue.put(new PairMersMap(k, tailMap, pairMersMap));
                }
            }
            mapsQueue.put(new PairMersMap(null));
        } catch (InterruptedException ex) {
            System.err.println(ex.getMessage());
        }

        //WAIT UNTIL DONE
        purgeMapsExecutorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            purgeMapsExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "PairMerSet purger interrupted exception!", TOOL_NAME);
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "PairMerSet purger execution exception!", TOOL_NAME);
            ex.printStackTrace();
            System.exit(1);
        } catch (TimeoutException ex) {
            ex.printStackTrace();
            Reporter.report("[ERROR]", "PairMerSet purger timeout exception!", TOOL_NAME);
        }
//        sizes = new int[mapChunks.size()];
//        for (int j = 0; j < mapChunks.size(); j++) {
////            System.err.println("chunk " + j + " size = " + mapChunks.get(j).size());
//            sizes[j] = mapChunks.get(j).size();
//            total += mapChunks.get(j).size();
//            min = sizes[j] < min ? sizes[j] : min;
//            max = sizes[j] > max ? sizes[j] : max;
////                        PairMersMap mapChunk = new PairMersMap(k, mapChunks.get(j), pairMersMap);
////                        mapsQueue.put(mapChunk);
//        }
//        if (sizes != null) {
//            Reporter.report("[INFO]", "Map split into " + sizes.length + " chunks, total elements=" + NumberFormat.getNumberInstance().format(total)
//                    + ", median=" + NumberFormat.getNumberInstance().format(CommonMaths.getMedian(sizes))
//                    + ", min=" + NumberFormat.getNumberInstance().format(min) + ", max=" + NumberFormat.getNumberInstance().format(max), TOOL_NAME);
//        }
        it = pairMerMaps.getkSizes().iterator();
        while (it.hasNext()) {
            Integer k = it.next();
            PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
//            System.err.println(k+"=k, |terminal|=" + pairMersMap.getTerminalPairMers().size());
            long size = pairMersMap.size();
            long terminals = pairMersMap.sizeTerminals();
            double ratio = (double) terminals / pairMersMap.getStoredSize();
            NumberFormat perc = NumberFormat.getPercentInstance();
            perc.setMaximumFractionDigits(4);
            Reporter.report("[INFO]", "Identified " + NumberFormat.getIntegerInstance().format(terminals) + " terminal PairMers (" + perc.format(ratio) + ") at k=" + k, TOOL_NAME);
            Reporter.report("[INFO]", "Found " + NumberFormat.getIntegerInstance().format(pairMersMap.getAmbiguous()) + " ambiguous extensions", TOOL_NAME);
            Reporter.report("[INFO]", "Finished purging map, k=" + k + ", n=" + NumberFormat.getIntegerInstance().format(size), TOOL_NAME);
//            Iterator<PairMer> iterator = pairMersMap.getTerminalPairMers().keySet().iterator();
//            while (iterator.hasNext()) {
//                PairMer next = iterator.next();
////                    System.err.println("terminal PM " + pairMersMap.get(next).getPairMerString(k, "_"));
//                if (!pairMersMap.contains(next)) {
//                    System.err.println("terminal PM absent from main map = " + pairMersMap.getTerminalPairMers().get(next).getPairMerString(k, "_"));
//                }
//            }
//            System.exit(1);
        }
    }

    private ConcurrentHashMap<Integer, PairMerToSeedMap> populateSeedMersMaps(ArrayList<Integer> kSizes) {

        ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeeds = new ConcurrentHashMap<>();

        //PREPARE EXECUTOR SERVICE
        int threads = kSizes.size();
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService endMersToSeedsExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //INIT THREADS
        for (Integer k : kSizes) {
            PairMerToSeedMapPopulator processor = new PairMerToSeedMapPopulator(seedSequences, kToSeeds, TOOL_NAME, k);
            futures.add(endMersToSeedsExecutorService.submit(processor));
        }

        //WAIT UNTIL DONE
        endMersToSeedsExecutorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            endMersToSeedsExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "SeedsProcessor interrupted exception!", TOOL_NAME);
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "SeedsProcessor execution exception!", TOOL_NAME);
//            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
            Reporter.report("[ERROR]", "SeedsProcessor timeout exception!", TOOL_NAME);
        }

//        for (Integer k : kSizes) {
//            Reporter.report("[INFO]", "Seed-mers map populated, k=" + k + ", n=" + seedSequences., TOOL_NAME);
//        }
        return kToSeeds;
    }

    private void purgeSeedMersFromPairMersMaps(ArrayList<Integer> kSizes, PairMerMaps pairMerMaps) {
        //PREPARE EXECUTOR SERVICE
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService seedMerPopulatorExecutorService = new ThreadPoolExecutor(threads + 1, threads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //INIT EMPTY SEED-PAIRMER-MAPS AND WRAP INPUT IN QUEUE
        PairMerMaps trimmedSeedPairMerMaps = new PairMerMaps(kSizes, TOOL_NAME);
        BlockingQueue<List<String>> seedsDummyQueue = new ArrayBlockingQueue<>(threads + 1);
        ArrayList<String> seedSequenceStrings = seedSequences.getSeedSequenceStrings();

        int size = seedSequenceStrings.size();
        int chunk = Math.max(1, (int) Math.ceil((double) size / threads));
//        System.err.println("Chunk="+chunk);
        try {
            for (int i = 0; i < size - chunk + 1; i += chunk) {
                seedsDummyQueue.put(seedSequenceStrings.subList(i, i + chunk));
//                System.err.println("ADDING elems "+i+" -> "+(i+chunk)+" of "+size);
            }
            seedsDummyQueue.put(new ArrayList<String>());
//            System.err.println("Finished adding elems ");
        } catch (InterruptedException ex) {

        }

        //INIT THREADS
        for (int i = 0; i < threads; i++) {
            PairMerMapPopulatorConsumer seedsPairMerMapPopulatorConsumer = new PairMerMapPopulatorConsumer(seedsDummyQueue,
                    trimmedSeedPairMerMaps, true, kSizes, MIN_KMER_FREQUENCY, true);
            futures.add(seedMerPopulatorExecutorService.submit(seedsPairMerMapPopulatorConsumer));
        }

        //WAIT UNTIL DONE
        seedMerPopulatorExecutorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            seedMerPopulatorExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "SeedsProcessor interrupted exception!", TOOL_NAME);
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "SeedsProcessor execution exception!", TOOL_NAME);
//            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
            Reporter.report("[ERROR]", "SeedsProcessor timeout exception!", TOOL_NAME);
        }

        for (Integer k : kSizes) {
//            Reporter.report("[INFO]", "Seed-mers map populated, k=" + k, TOOL_NAME);

            Iterator<PairMer> it = trimmedSeedPairMerMaps.getPairMersMap(k).iterator();
            PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);

            long removedCount = 0L;
            while (it.hasNext()) {
                if (pairMersMap.remove(it.next())) {
                    removedCount++;
                }
            }
            Reporter.report("[INFO]", "Seed-based purging removed additional " + NumberFormat.getIntegerInstance().format(removedCount) + " for k=" + k, TOOL_NAME);
        }
    }

    private void extendSeedsParallel(ArrayList<Integer> kSizes, PairMerMaps pairMerMaps,
            ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers, OptSet optSet) {
        //PREPARE EXECUTOR SERVICE
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService seedExtenderExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        BlockingQueue<PairMersMap> queue = new ArrayBlockingQueue<>(kSizes.size() + 1);
        try {
            for (int k : kSizes) {
                queue.put(pairMerMaps.getPairMersMap(k));
            }
            queue.put(new PairMersMap(null));
        } catch (InterruptedException ex) {
        }

        byte threadId = Byte.MIN_VALUE;
        //INIT THREADS
        for (int i = 0; i < threads; i++) {
//            PairMersSeedExtenderConsumer = new PairMersSeedExtenderConsumer
            PairMersSeedExtenderConsumer extenderConsumer = new PairMersSeedExtenderConsumer(queue, kToSeedMers,
                    (String) optSet.getOpt("debug-file").getValueOrDefault(), TOOL_NAME, threadId++);
            futures.add(seedExtenderExecutorService.submit(extenderConsumer));
        }

        //WAIT UNTIL DONE
        seedExtenderExecutorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            seedExtenderExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "SeedsExtender interrupted exception!", TOOL_NAME);
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "SeedsExtender execution exception!", TOOL_NAME);
//            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
            Reporter.report("[ERROR]", "SeedsExtender timeout exception!", TOOL_NAME);
        }

    }

    private int setKmerLength(int k) {
        KMER_LENGTH_MIN = k;
        KMER_LENGTH_MAX = k;
        KMER_LENGTH_STEP = 1;
        return k;
    }

//
//    private void printHelp() {
//
//        System.out.println();
//        System.out.println("Given a set of k-mers, performs their unmabiguous extension.");
//        System.out.println("Input may be passed through stdin. Output is printed to stdout. ");
//        System.out.println("Progress information, warnings errors etc. are printed to stderr, usage:");
//        System.out.println("java -jar " + this.TOOL_NAME + ".jar [options] [INPUT_FILENAME(s)] ");
//        System.out.println(" -h, --help                     : Print this help screen and exit");
//        System.out.println(" -k, --kmer-length <int>        : Required only if input other than a list of k-mers");
//        System.out.println(" -t, --threads <int>            : Number of threads, defaults to " + MAX_THREADS);
//        System.out.println("                                : ");
//        System.out.println(" -F, --fasta-out                : Output each k-mer as a separate FASTA record instead of just listing extended nucleotide sequences");
//        System.out.println(" -N, --name-prefix <string>     : If outputing FASTA, prefix each record identifier with this <string>");
//        System.out.println("                                : ");
//        System.out.println(" -O, --out-redirect <string>    : Redirect stdout to this file");
//        System.out.println(" -E, --err-redirect <string>    : Redirect stderr to this file");
//        System.out.println(" -S, --stats-file <string>      : Write extension stats to this file");
//        System.out.println(" -D, --debug-file <string>      : Write unkosher extensions details to this file");
//        System.out.println("                                : ");
////        System.out.println(" Multi-pass hashing options, incompatible with stdin:");
////        System.out.println(" -M, --multipass-compress <int> : Reduce memory requirements with <int> hashing iterations");
////        System.out.println(" -H, --hash-array-size    <int> : At each iteration, array of <int> size will be added, defualts to " + NumberFormat.getNumberInstance().format(HASH_ARRAY_SIZE) + "*");
////        System.out.println("                                : ");
////
////        String s = "* - Note that each hash array will consume 4 Bytes per element, so the largest possible (on most VMs) "
////                + "2,147,483,644 element array will consume over 8GB of the allocated memory. ";
////        System.out.println(Reporter.wrapString(s, 145));
//        String s = "Currently k-mer frequency is not taken into consideration, so use of a dedicated k-mer counting program, "
//            + "such as KMC or Jellyfish is recommended. It is best to exclude low frequency k-mers before passing "
//            + "the list of k-mers to KmerExtender. For smaller jobs FASTA or FASTQ input may suffice.";
//        System.out.println(Reporter.wrapString(s, 145));
//        System.out.println();
//    }
}
