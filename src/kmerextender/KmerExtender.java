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
package kmerextender;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import shared.Reporter;
import shared.InputReaderProducer;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.print.Collation;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerExtender {

//    private PairMersMap pairMersMap = new PairMersMap(); //TODO user to give expected num of elems?
    private final ArrayList<String> inputFileNamesList = new ArrayList<>();
//    private Integer KMER_LENGTH;
    private Integer KMER_LENGTH_MIN;
    private Integer KMER_LENGTH_MAX;
    private Integer KMER_LENGTH_STEP;
    private SeedSequences seedSequences;
    private Integer MIN_KMER_FREQUENCY;
    private int MAX_THREADS;
    private int INPUT_BUFFER_SIZE;
    private int INPUT_QUEUE;

    private boolean OUTPUT_FASTA = false;
    private String NAME_PREFIX = "";
    private String DEBUG_FILE;
    private String STATS_FILE;
//    private Integer HASH_ARRAY_SIZE = Integer.MAX_VALUE - 3;
//    private Integer MULTIPASS_COMPRESS = null;

    private final String TOOL_NAME;
    private final int HELP_WIDTH = 170;

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
        runKmerExtender();
//        }
    }

    private void readArgValues(OptSet optSet) {
        INPUT_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        INPUT_QUEUE = (int) optSet.getOpt("Q").getValueOrDefault();

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
            seedSequences = new SeedSequences((String) optSet.getOpt("seed-file").getValueOrDefault());
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
        String outRedirect;
        String errRedirect;
        if ((outRedirect = (String) optSet.getOpt("o").getValueOrDefault()) != null) {
            try {
                File file = new File(outRedirect);
                PrintStream printStream;
                printStream = new PrintStream(new FileOutputStream(file));
                System.setOut(printStream);
            } catch (FileNotFoundException ex) {
                Reporter.report("[ERROR]", "Failed redirecting stdout to " + outRedirect, TOOL_NAME);
            }
        }
        if ((errRedirect = (String) optSet.getOpt("e").getValueOrDefault()) != null) {
            try {
                File file = new File(errRedirect);
                PrintStream printStream;
                printStream = new PrintStream(new FileOutputStream(file));
                System.setErr(printStream);
            } catch (FileNotFoundException ex) {
                Reporter.report("[ERROR]", "Failed redirecting stderr to " + errRedirect, TOOL_NAME);
            }
        }

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
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input settings - general extender]");
        optSet.addOpt(new Opt('k', "k-mer-length", "Required only if input other than a list of k-mers", 1).setMinValue(4).setMaxValue(2048));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of records (k-mers or FASTQ reads or pairs depending on input) "
            + "passed to in-queue", 1024, 128, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
            2, 1, 256));
//        optSet.addOpt(new Opt('m', "min-frequency", "....", 1,1,Integer.MAX_VALUE));
        int footId = 1;
        String foot = "Warning! Exploring a large range of k values for a significant input [k-mers/FAST[A|Q]] "
            + "will make the memory requirements explode, "
            + "as it is done in parallel to avoid excessive I/O and to preserve stdin handling";
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Variable k-mer size for longest extension of a seed]");
        optSet.addOpt(new Opt('S', "seed-file", "Fasta file containing \"seeds\" to be extended", 1));
        optSet.addOpt(new Opt(null, "k-mer-min", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-max", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-step", "", 1).setMinValue(1).setDefaultValue(2).addFootnote(footId, foot));
//        optSet.addOpt(new Opt(null, "k-values", "Explicitly set k values to be explored", 2).setMinValue(1).setDefaultValue(2).addFootnote(footId, foot));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('t', "threads", "Number threads for populating a set of implicitly connected k-mers",
            1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('f', "fasta-out", "Output each k-mer as a separate FASTA record instead of just listing extended nucleotide sequences"));
        optSet.addOpt(new Opt('p', "fasta-id-prefix", "Prefix each FASTA identifier with <arg> ", 1));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
        optSet.addOpt(new Opt('s', "stats-file", "Write extension stats to this file", 1));
        optSet.addOpt(new Opt('d', "debug-file", "Write unkosher extensions details to this file", 1));

        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    /**
     * Wrapper method executing individual steps
     */
    private void runKmerExtender() {
        Reporter.report("[INFO]", "Initialized, will use " + MAX_THREADS + " thread(s) to populate map ", TOOL_NAME);

        ArrayList<Integer> kSizes = new ArrayList<>();
        if (KMER_LENGTH_MIN != null) {
            for (int k = KMER_LENGTH_MIN; k < KMER_LENGTH_MAX + 1; k += KMER_LENGTH_STEP) {
                kSizes.add(k);
            }
        }
        Collections.sort(kSizes, Collections.reverseOrder()); 

        //READ k-mers AND POPULATE A MAP FOR EACH SIZE OF k
        PairMerMaps pairMerMaps = new PairMerMaps(kSizes);
        readKmersAndPopulatePairMersMaps(pairMerMaps);
        for (Integer k : kSizes) {
            PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
            Reporter.report("[INFO]", "Finished populating map, k=" + k + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);
        }

        //RELEASE SOME MEMORY
        Reporter.report("[INFO]", "Now purging to release some memory...", TOOL_NAME);
        purgePopulatedPairMersMaps(pairMerMaps);

        ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers = null;
        if (seedSequences != null) {
            //PURGE NON-TERMINAL SeedPairMers FROM PairMerMaps 
            purgeSeedMersFromPaierMersMaps(kSizes, pairMerMaps);
            //CREATE PairMers REPRESENTING SEED-ENDS 
            kToSeedMers = populateSeedMersMaps(kSizes);
        }

        //EXTEND - CAN BE PARALLELIZED IF NO INTENTION TO GENERATE CROSS-k-EXTENSIONS 
        if (seedSequences == null) {
            for (Integer k : kSizes) {
                PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);
                PairMersExtender pairMersExtender = new PairMersExtender(DEBUG_FILE, STATS_FILE, TOOL_NAME);
                pairMersExtender.matchAndExtendKmers(k, pairMersMap, OUTPUT_FASTA, NAME_PREFIX);
                Reporter.report("[INFO]", "Finished extending for k=" + k, TOOL_NAME);

            }
//            pairMersExtender.matchAndExtendSeeds(k, pairMersMap, kToSeedMers.get(k));
        } else {
            extendSeedsParallel(kSizes, pairMerMaps, kToSeedMers);

        }

        //SELECT LONGEST EXTENSION
        if (seedSequences != null) {
            for (SeedSequence seed : seedSequences.getSeedSequences()) {
//            Reporter.report("[INFO]", "Longest extension of the provided seed is from " + seed.getSequenceString().length() + " to " + longestSeedAfterExtension + " at k = " + kBest, TOOL_NAME);
//                Map.Entry<Integer, String> longest = seed.getLongestExtended();
                Map.Entry<Integer, SeedExtensionsPair> longestExtensionLeft = seed.getLongestExtensionLeft();
                Map.Entry<Integer, SeedExtensionsPair> longestExtensionRight = seed.getLongestExtensionRight();
//                StringBuilder output = new StringBuilder();
                StringBuilder outputLR = new StringBuilder();
                String extensionLeft = longestExtensionLeft.getValue().getExtensionLeft();
                String extensionRight = longestExtensionRight.getValue().getExtensionRight();
                if (OUTPUT_FASTA) {
//                    if (!NAME_PREFIX.isEmpty()) {
//                        output.append(">").append(NAME_PREFIX).append(" ");
//                        output.append(longest.getValue().length()).append(System.lineSeparator());
//                    } else if (longest.getKey() == 0) { //NOEXTENSION BEYOND INPUT                        
//                        output.append(">").append(seed.getId());
//                        output.append(" not_extended=").append(longest.getValue().length());
//                        output.append(System.lineSeparator());
//                    } else {
//                        output.append(">").append(seed.getId());
//                        output.append(" extended from=").append(seed.getSequenceString().length());
//                        output.append(" to=").append(longest.getValue().length()).append(" at k=").append(longest.getKey());
//                        output.append(System.lineSeparator());
//                    }

                    if (!NAME_PREFIX.isEmpty()) {
                        outputLR.append(">").append(NAME_PREFIX);
                    } else {
                        outputLR.append(">").append(seed.getId());
                    }
                    outputLR.append(" L=").append(extensionLeft.length()).append(" at k=").append(longestExtensionLeft.getKey());
                    outputLR.append(" R=").append(extensionRight.length()).append(" at k=").append(longestExtensionRight.getKey());
                    outputLR.append(" extended=");
                    outputLR.append(extensionLeft.length() + seed.getSequenceString().length() + extensionRight.length());
                    outputLR.append(System.lineSeparator());

                }
//                output.append(longest.getValue());
//                System.out.println(output);

                outputLR.append(extensionLeft).append(seed.getSequenceString()).append(extensionRight);
                System.out.println(outputLR);

            }
        }
        Reporter.report("[INFO]", "Finished extending k-mers", TOOL_NAME);
    }

    /**
     * Producer-consumer multi-threading to: - read the input (k-mers, FASTA, etc) - generate PairMer objects - populate
     * Map(s)
     */
    private void readKmersAndPopulatePairMersMaps(PairMerMaps pairMerMaps) {
        BlockingQueue inputQueue = new ArrayBlockingQueue(INPUT_QUEUE);

        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = MAX_THREADS;
            Reporter.report("[INFO]", "Allocated " + threads + " thread(s) to map populating", TOOL_NAME);
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateExecutor = new ThreadPoolExecutor(threads+1, threads+1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

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
            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
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
        //PREPARE QUEUE OF MAPS TO BE PICKED-UP BY THREADS
        ArrayBlockingQueue<PairMersMap> mapsQueue = new ArrayBlockingQueue(pairMerMaps.size() + 1);
        Iterator<Integer> it = pairMerMaps.getkSizes().iterator();
        try {
            while (it.hasNext()) {
                mapsQueue.put(pairMerMaps.getPairMersMap(it.next()));
            }
            mapsQueue.put(new PairMersMap(null));
        } catch (InterruptedException ex) {
            System.err.println(ex.getMessage());
        }
        //PREPARE EXECUTOR SERVICE
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService purgeMapsExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        //INIT THREADS
        for (int i = 0; i < threads; i++) {
            futures.add(purgeMapsExecutorService.submit(new PairMerMapPurger(mapsQueue, TOOL_NAME)));
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
//            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
            Reporter.report("[ERROR]", "PairMerSet purger timeout exception!", TOOL_NAME);
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

    private void purgeSeedMersFromPaierMersMaps(ArrayList<Integer> kSizes, PairMerMaps pairMerMaps) {
        //PREPARE EXECUTOR SERVICE        
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService seedMerPopulatorExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //INIT EMPTY SEED-PAIRMER-MAPS AND WRAP INPUT IN QUEUE
        PairMerMaps trimmedSeedPairMerMaps = new PairMerMaps(kSizes);
        BlockingQueue<List<String>> seedsDummyQueue = new ArrayBlockingQueue<>(threads + 1);
        ArrayList<String> seedSequenceStrings = seedSequences.getSeedSequenceStrings();

        int size = seedSequenceStrings.size();
        int chunk = Math.max(1, size / threads);
        try {
            for (int i = 0; i < size - chunk + 1; i += chunk) {
                seedsDummyQueue.put(seedSequenceStrings.subList(i, i + chunk));
//                System.err.println("ADDING elems "+i+" -> "+(i+chunk)+" of "+size);
            }
            seedsDummyQueue.put(new ArrayList<String>());
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

            Iterator<PairMer> it = trimmedSeedPairMerMaps.getPairMersMap(k).getPairMersSkipListMap().keySet().iterator();
            PairMersMap pairMersMap = pairMerMaps.getPairMersMap(k);

            long removedCount = 0L;
            while (it.hasNext()) {
                if (pairMersMap.getPairMersSkipListMap().remove(it.next()) != null) {
                    removedCount++;
                }
            }
            Reporter.report("[INFO]", "Seed-based purging removed additional " + NumberFormat.getIntegerInstance().format(removedCount), TOOL_NAME);
        }

//            PairMersMap pairMersMap = kToMap.get(k);
//            PairMerToSeedMap pairMerToSeedMap = new PairMerToSeedMap(seedSequences, k);
//            kToSeeds.put(actualK, pairMerToSeedMap);
//            BlockingQueue<ArrayList<String>> seedsDummyQueue = new ArrayBlockingQueue<>(2);
//            PairMersMap trimmedSeedPairMersMap = new PairMersMap();
//            try {
//                seedsDummyQueue.put(seedSequences.getSeedSequenceStrings());
//                seedsDummyQueue.put(new ArrayList<String>());
//            } catch (InterruptedException ex) {
//            }
//            PairMerMapPopulatorConsumer seedsPairMerMapPopulatorConsumer = new PairMerMapPopulatorConsumer(seedsDummyQueue,
//                trimmedSeedPairMersMap, true, k, MIN_KMER_FREQUENCY, true);
//            seedsPairMerMapPopulatorConsumer.run(); //TODO move to a background thread earlier  if needed
//            Reporter.report("[INFO]", "Seed-mers map populated, k=" + actualK, TOOL_NAME);
//            Iterator<PairMer> it = trimmedSeedPairMersMap.getPairMersSkipListMap().keySet().iterator();
//            long removedCount = 0L;
//            while (it.hasNext()) {
//                if (pairMersMap.getPairMersSkipListMap().remove(it.next()) != null) {
//                    removedCount++;
//                }
//            }
//            Reporter.report("[INFO]", "Seed-based purging removed additional " + NumberFormat.getIntegerInstance().format(removedCount), TOOL_NAME);
//        }
//        //PREPARE QUEUE OF MAPS TO BE PICKED-UP BY THREADS
//        ArrayBlockingQueue<PairMersMap> mapsQueue = new ArrayBlockingQueue(kToMap.size()+1);
//        Iterator<PairMersMap> it = kToMap.values().iterator();
//        try {
//            while (it.hasNext()) {
//                mapsQueue.put(it.next());
//            }
//            mapsQueue.put(new PairMersMap());
//        } catch (InterruptedException ex) {
//            System.err.println(ex.getMessage());
//        }
//        //PREPARE EXECUTOR SERVICE
//        int threads = MAX_THREADS;
//        ArrayList<Future<?>> futures = new ArrayList<>(threads);
//        final ExecutorService purgeMapsExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
//        //INIT THREADS
//        for (int i = 0; i < threads; i++) {
//            futures.add(purgeMapsExecutorService.submit(new PairMerMapPurger(mapsQueue)));
//        }
//        //WAIT UNTIL DONE
//        purgeMapsExecutorService.shutdown();
//        try {
//            for (Future<?> f : futures) {
//                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            }
//            purgeMapsExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            Reporter.report("[ERROR]", "PairMerSet purger interrupted exception!", TOOL_NAME);
//        } catch (ExecutionException ex) {
//            Reporter.report("[ERROR]", "PairMerSet purger execution exception!", TOOL_NAME);
////            ex.printStackTrace();
//        } catch (TimeoutException ex) {
//            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
//            Reporter.report("[ERROR]", "PairMerSet purger timeout exception!", TOOL_NAME);
//        }
    }

    private void extendSeedsParallel(ArrayList<Integer> kSizes, PairMerMaps pairMerMaps, ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers) {
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

        //INIT THREADS
        for (int i = 0; i < threads; i++) {
//            PairMersSeedExtenderConsumer = new PairMersSeedExtenderConsumer
            PairMersSeedExtenderConsumer extenderConsumer = new PairMersSeedExtenderConsumer(queue, kToSeedMers, DEBUG_FILE, TOOL_NAME);
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
