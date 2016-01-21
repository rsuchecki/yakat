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
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.Sequence;

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
    private SeedSequences SEED_SEQUENCES;
    private Integer MIN_KMER_FREQUENCY;
    private int MAX_THREADS;
    private boolean OUTPUT_FASTA = false;
    private String NAME_PREFIX = "";
    private String DEBUG_FILE;
    private String STATS_FILE;
//    private Integer HASH_ARRAY_SIZE = Integer.MAX_VALUE - 3;
//    private Integer MULTIPASS_COMPRESS = null;

    private boolean RUN_SOME_WILD_AND_WONDERFUL_STUFF = false;
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
        if (RUN_SOME_WILD_AND_WONDERFUL_STUFF) {
//            new kmerextender.ideas.Alternative(MAX_THREADS, inputFileNamesList, KMER_LENGTH);
        } else {
            runKmerExtender();
        }
    }

    private void readArgValues(OptSet optSet) {
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
            SEED_SEQUENCES = new SeedSequences((String) optSet.getOpt("seed-file").getValueOrDefault());
            if (SEED_SEQUENCES.getSeedSequences().isEmpty()) {
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
                + "as it is done in parallel to avoid excessive I/O and preserve stdin handling";
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Variable k-mer size for longest extension of a single seed]");
        optSet.addOpt(new Opt('S', "seed-file", "Fasta file containing a single entry (\"a seed\") to be extended", 1));
        optSet.addOpt(new Opt(null, "k-mer-min", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-max", "", 1).setMinValue(3).setMaxValue(2048).addFootnote(footId, foot));
        optSet.addOpt(new Opt(null, "k-mer-step", "", 1).setMinValue(1).setDefaultValue(2).addFootnote(footId, foot));
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

    private void runKmerExtender() {
        Reporter.report("[INFO]", "Initialized, will use " + MAX_THREADS + " thread(s) to populate map ", TOOL_NAME);

        HashMap<Integer, PairMersMap> kSizeToPairMersMap = new HashMap<>();
        ArrayList<Integer> kSizes = new ArrayList<>();
        if (KMER_LENGTH_MIN != null) {
            for (int k = KMER_LENGTH_MIN; k < KMER_LENGTH_MAX + 1; k += KMER_LENGTH_STEP) {
                kSizeToPairMersMap.put(k, new PairMersMap());
                kSizes.add(k);
            }
        } else {
            kSizes.add(0); // no k-size given, will be taken from input k-mers or the process will terminate if FASTA/FASTQ input
            kSizeToPairMersMap.put(0, new PairMersMap());
        }

        //READ k-mers AND POPULATE A MAP FOR EACH SIZE OF k
        int actualK = readKmersAndPopulatePairMersMaps(kSizeToPairMersMap, kSizes);

        for (Integer k : kSizes) {
            actualK = k > 0 ? k : actualK;            
            PairMersMap pairMersMap = kSizeToPairMersMap.get(k);
            Reporter.report("[INFO]", "Finished populating map, k=" + actualK + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);
        }
        
        int threads = MAX_THREADS;
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        final ExecutorService purgeMapAndExtendExecutorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        
        for (Integer k : kSizes) {
            actualK = k > 0 ? k : actualK;            
            PairMersMap pairMersMap = kSizeToPairMersMap.get(k);
            long purged = pairMersMap.purge(k); //indicating large number of kmers overall        
            if (purged > 100000) {
                for (int i = 0; i < 5; i++) {
                    System.gc();
                    try {
                        Thread.sleep(500);
                    } catch (InterruptedException ex) {
                    }
                }
            }
            Reporter.report("[INFO]", "Finished purging map, k= " + actualK + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);

            //Seed-processing ================================== BEGIN ============================
            if (SEED_SEQUENCES != null) {
                BlockingQueue<ArrayList<String>> seedsDummyQueue = new ArrayBlockingQueue<>(2);
                ArrayList<String> seedStrings = new ArrayList<>(1);
                for (SeedSequence seed : SEED_SEQUENCES.getSeedSequences()) {
                    String trimmedSeed = seed.getSequenceString().substring(1, seed.getSequenceString().length() - 1); //hack to prevent end-pairmers from being thrown out 
                    seedStrings.add(trimmedSeed);
                }
                PairMersMap trimmedSeedPairMersMap = new PairMersMap();
                try {
                    seedsDummyQueue.put(seedStrings);
                    seedsDummyQueue.put(new ArrayList<String>());
                } catch (InterruptedException ex) {
                }
                PairMerMapPopulatorConsumer seedsPairMerMapPopulatorConsumer = new PairMerMapPopulatorConsumer(seedsDummyQueue, trimmedSeedPairMersMap,
                        true, k, MIN_KMER_FREQUENCY);
                seedsPairMerMapPopulatorConsumer.run(); //TODO move to a background thread earlier  if needed
                Iterator<PairMer> it = trimmedSeedPairMersMap.getPairMersSkipListMap().keySet().iterator();
                long removedCount = 0L;
                while (it.hasNext()) {
                    if (pairMersMap.getPairMersSkipListMap().remove(it.next()) != null) {
                        removedCount++;
                    }
                }
                Reporter.report("[INFO]", "Seed-based purging removed additional " + NumberFormat.getIntegerInstance().format(removedCount), TOOL_NAME);
            }
            //Seed-processing ================================== END ============================

            
            
            PairMersExtender pairMersExtender = new PairMersExtender(DEBUG_FILE, STATS_FILE, TOOL_NAME);
            pairMersExtender.matchAndExtendKmers(actualK, pairMersMap, OUTPUT_FASTA, NAME_PREFIX, MAX_THREADS, SEED_SEQUENCES);
            
//            futures.add(purgeMapAndExtendExecutorService.submit(consumer));
            
        }
        
//        purgeMapAndExtendExecutorService.shutdown();
//        try {
//            for (Future<?> f : futures) {
//                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            }
//            purgeMapAndExtendExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            Reporter.report("[ERROR]", "purgeMapAndExtendExecutorService interrupted exception!", TOOL_NAME);
//        } catch (ExecutionException ex) {
//            Reporter.report("[ERROR]", "purgeMapAndExtendExecutorService execution exception!", TOOL_NAME);
//            ex.printStackTrace();
//
//        } catch (TimeoutException ex) {
//            Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
//            Reporter.report("[ERROR]", "purgeMapAndExtendExecutorService timeout exception!", TOOL_NAME);
//        }
        
        
        
        if (SEED_SEQUENCES != null) {
            for (SeedSequence seed : SEED_SEQUENCES.getSeedSequences()) {
//            Reporter.report("[INFO]", "Longest extension of the provided seed is from " + seed.getSequenceString().length() + " to " + longestSeedAfterExtension + " at k = " + kBest, TOOL_NAME);
                Map.Entry<Integer, String> longest = seed.getLongestExtended();
                if (OUTPUT_FASTA) {
                    if (!NAME_PREFIX.isEmpty()) {
                        System.out.println(">" + NAME_PREFIX + " " + longest.getValue().length());
                    } else if (longest.getKey() == 0) { //NOEXTENSION BEYOND INPUT
                        System.out.println(">" + seed.getId() + " not_extended=" + longest.getValue().length());
                    } else {
                        System.out.println(">" + seed.getId() + " extended from=" + seed.getSequenceString().length()
                                + " to=" + longest.getValue().length() + " at k=" + longest.getKey());
                    }

                }
                System.out.println(longest.getValue());

            }
        }
        Reporter.report("[INFO]", "Finished extending k-mers", TOOL_NAME);
    }

    /**
     * Producer-consumer multi-threading to read input (k-mers, FASTA, etc)
     * generate PairMer objects and populate a Map
     */
    private int readKmersAndPopulatePairMersMaps(HashMap<Integer, PairMersMap> kSizeToPairMersMap,
            ArrayList<Integer> kValues) {
        int kReturn = 0;

        HashMap<Integer, BlockingQueue> kSizeToinputQueue = new HashMap<>();

        //POSSIBLY single 0 value when k is to be picked-up from input (if list of k-mers given as input)
        for (Integer kValue : kValues) {
            kSizeToinputQueue.put(kValue, new ArrayBlockingQueue(1024));
        }

        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = MAX_THREADS;
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateDbExecutor = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            //SPAWN INPUT READING THREAD
            InputReaderProducer inputReaderProducer = new InputReaderProducer(kSizeToinputQueue, kValues, inputFileNamesList, TOOL_NAME);

            Future<?> future = readAndPopulateDbExecutor.submit(inputReaderProducer);
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
            Integer kSizeFromInput = null;
            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
                splitInputSequenceIntoKmers = false;
                kSizeFromInput = setKmerLength(inputReaderProducer.getKmerLength());
            }

            //SPAWN THREADS TO POPULATE MAP
            for (Integer kValue : kValues) {
                BlockingQueue queue = kSizeToinputQueue.get(kValue);
                PairMersMap pairMersMap = kSizeToPairMersMap.get(kValue);
                if (kValue == 0) {   //single dummy-entry
                    kValue = kSizeFromInput;
                    kReturn = kSizeFromInput;
                }
                threads = Math.max(1, threads / kValues.size()); //REDUCE NUM_THREADS BASED ON NUMER OF k-values To BE INVESTIGATED
                for (int i = 0; i < threads; i++) {
                    PairMerMapPopulatorConsumer consumer = new PairMerMapPopulatorConsumer(queue, pairMersMap, splitInputSequenceIntoKmers, kValue, MIN_KMER_FREQUENCY);
//                    PairMerMapPopulatorConsumer consumer = new PairMerMapPopulatorConsumer(inputQueue, pairMersMap, splitInputSequenceIntoKmers, KMER_LENGTH, MIN_KMER_FREQUENCY);
                    futures.add(readAndPopulateDbExecutor.submit(consumer));
                }
            }
            readAndPopulateDbExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                readAndPopulateDbExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
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

        return kReturn;

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
