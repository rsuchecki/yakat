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
        }

        //READ k-mers AND POPULATE A MAP FOR EACH SIZE OF k 
        readKmersAndPopulatePairMersMaps(kSizeToPairMersMap, kSizes);

        for (Integer k : kSizes) {
            PairMersMap pairMersMap = kSizeToPairMersMap.get(k);
            Reporter.report("[INFO]", "Finished populating map, k=" + k + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);
        }
        for (Integer k : kSizes) {
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
            Reporter.report("[INFO]", "Finished purging map, k= " + k + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);

            //Seed-processing ================================== BEGIN ============================
            if (SEED_SEQUENCES != null) {
                BlockingQueue<ArrayList<String>> seedsDummyQueue = new ArrayBlockingQueue<>(2);
                ArrayList<String> seedStrings = new ArrayList<>(1);
                for(SeedSequence seed: SEED_SEQUENCES.getSeedSequences()) {
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
            pairMersExtender.matchAndExtendKmers(k, pairMersMap, OUTPUT_FASTA, NAME_PREFIX, MAX_THREADS, SEED_SEQUENCES);
//            if (seedPossiblyExtended != null) {
//                if (longestSeedAfterExtension < seedPossiblyExtended.length()) {
//                    longestSeedAfterExtension = seedPossiblyExtended.length();
////                    kBest = k;
//                    seedToOutput = seedPossiblyExtended;
//                }
////                System.out.println(seedPossiblyExtended.length() + " @ k = " + k);
//            }
        }
        for(SeedSequence seed: SEED_SEQUENCES.getSeedSequences()) {        
//            Reporter.report("[INFO]", "Longest extension of the provided seed is from " + seed.getSequenceString().length() + " to " + longestSeedAfterExtension + " at k = " + kBest, TOOL_NAME);
            Map.Entry<Integer, String> longest = seed.getLongestExtended();            
            if (OUTPUT_FASTA) {
                if (!NAME_PREFIX.isEmpty()) {
                    System.out.println(">" + NAME_PREFIX + " " + longest.getValue().length());
                } else if(longest.getKey() == 0){ //NOEXTENSION BEYOND INPUT
                    System.out.println(">" + seed.getId() + " not extended " + longest.getValue().length());
                } else {
                    System.out.println(">" + seed.getId() + " extended from " + seed.getSequenceString().length() + 
                            " to " + longest.getValue().length() + " at k = " + longest.getKey());
                }

//                } else {
//                    System.out.println(">InputSeedID " + inputSeedLen);
//                }
//                System.out.println(seed);            
            }
            System.out.println(longest.getValue());

        }
        Reporter.report("[INFO]", "Finished extending k-mers", TOOL_NAME);
    }

    /**
     * Producer-consumer multi-threading to read input (k-mers, FASTA, etc)
     * generate PairMer objects and populate a Map
     */
    private HashMap<Integer, PairMersMap> readKmersAndPopulatePairMersMaps(HashMap<Integer, PairMersMap> kSizeToPairMersMap,
            ArrayList<Integer> kValues) {

        HashMap<Integer, BlockingQueue> kSizeToinputQueue = new HashMap<>();

        //POSSIBLY single 0 value when k is to be picked-up from input (if list of k-mers given)
        for (Integer kValue : kValues) {
            kSizeToinputQueue.put(kValue, new ArrayBlockingQueue(1024));
        }

        try {
//            Integer k;
            //READ INPUT AND POPULATE PairMers MAP
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);
            int threads = MAX_THREADS;
//            BlockingQueue inputQueue = new ArrayBlockingQueue(65536);
//            BlockingQueue inputQueue = new ArrayBlockingQueue(1024);
//            boolean stranded = false;
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateDbExecutor = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            //SPAWN INPUT READING THREAD
            InputReaderProducer inputReaderProducer = new InputReaderProducer(kSizeToinputQueue, kValues, inputFileNamesList, TOOL_NAME);
//                inputReaderProducer = new InputReaderProducer(inputQueue, inputFileNamesList, KMER_LENGTH, pairMersMap, TOOL_NAME);

            Future<?> future = readAndPopulateDbExecutor.submit(inputReaderProducer);
            futures.add(future);
            boolean splitInputSequenceIntoKmers = true;
//            if(KMER_LENGTH == null) {
//               splitInputSequenceIntoKmers = false;
//            }

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
                if (kValue == 0) {   //single dummy-entry
                    kValue = kSizeFromInput;
                }
                threads = Math.max(1, threads / kValues.size()); //REDUCE NUM_THREADS BASED ON NUMER OF k-values To BE INVESTIGATED
                for (int i = 0; i < threads; i++) {
                    PairMerMapPopulatorConsumer consumer = new PairMerMapPopulatorConsumer(queue, kSizeToPairMersMap.get(kValue), splitInputSequenceIntoKmers, kValue, MIN_KMER_FREQUENCY);
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
                ex.getMessage();

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

        return kSizeToPairMersMap;

    }

    private int setKmerLength(int k) {
        KMER_LENGTH_MIN = k;
        KMER_LENGTH_MAX = k;
        KMER_LENGTH_STEP = 1;
        return k;
    }

//    private int getKmerLength() {
//        if (KMER_LENGTH_MIN == KMER_LENGTH_MAX) {
//            return KMER_LENGTH_MIN;
//        } else {            
//            return -1; 
//        }
//    }
//    /**
//     * Supposedly POSIX-compliant CLI args processor
//     *
//     * @param args
//     */
//    private void processArgs(String[] args) {
//        ArrayList<String> allArgs = new ArrayList<>();
//        //SPLIT-UP POSIX STYLE (-abcd) IF ANY
//        for (String arg : args) {
//            if (arg.startsWith("--")) {
//                allArgs.add(arg);
//            } else if (arg.startsWith("-")) {
//                if (arg.length() > 2) {
//                    char[] chars = arg.toCharArray();
//                    StringBuilder sb = new StringBuilder();
//                    for (int i = 1; i < chars.length; i++) {
//                        char c = chars[i];
////                        System.out.println(c + " ASCII value " + (int) c);
//                        if ((int) c >= 48 && (int) c <= 57) {
//                            sb.append(c);
//                        } else {
//                            allArgs.add("-" + c);
//                        }
//                    }
//                    String toString = sb.toString();
//                    if (!toString.isEmpty()) {
//                        allArgs.add(toString);
//                    }
//
//                } else {
//                    allArgs.add(arg);
//                }
//            } else if (!arg.isEmpty()) { //we get en empty oone after tool/module name has been removed by outer wrapper
//                allArgs.add(arg);
//            }
//        }
//        //Now process args
//        Iterator<String> it = allArgs.iterator();
//
//        while (it.hasNext()) {
//            String a = it.next();
//            //
//            if (a.equals("-h") || a.equals("--help")) {
//                printHelp();
//                System.exit(0);
//            } else if (a.equals("-k") || a.equalsIgnoreCase("--kmer-length")) {
//                if (it.hasNext()) {
//                    String optionValue = it.next();
//                    try {
//                        KMER_LENGTH = Integer.parseInt(optionValue);
//                        if (KMER_LENGTH < 4) {
//                            System.err.println("Fatal error! \"k\" must be set to at least 4, currently set to: " + KMER_LENGTH);
//                            System.exit(1);
//                        }
//                    } catch (NumberFormatException e) {
//                        System.err.println("Fatal error! The k-mer lenght must be given as an integer, offending argument: \"" + optionValue + "\"");
//                    }
//                }
//            } else if (a.equals("-t") || a.equalsIgnoreCase("--threads")) {
//                if (it.hasNext()) {
//                    String optionValue = it.next();
//                    try {
//                        MAX_THREADS = Integer.parseInt(optionValue);
//                        if (MAX_THREADS < 1) {
//                            System.err.println("Fatal error! \"threads\" must be set to a positive integer, currently set to: " + MAX_THREADS);
//                            System.exit(1);
//                        } else if (MAX_THREADS > Runtime.getRuntime().availableProcessors()) {
//                            MAX_THREADS = Runtime.getRuntime().availableProcessors();
//                            System.err.println("Warning, number of threads reset to the maximum available (" + MAX_THREADS + ") CPU cores, offending argument: \"" + optionValue + "\"");
//                        }
//                    } catch (NumberFormatException e) {
//                        System.err.println("Fatal error! Number of threads must be given as an integer, offending argument: \"" + optionValue + "\"");
//                    }
//                }
//            } else if (a.equals("-D") || a.equalsIgnoreCase("--debug-file")) {
//                if (it.hasNext()) {
//                    DEBUG_FILE = it.next();
//                    Reporter.writeToFile(DEBUG_FILE, Reporter.formatReport("[DEBUG]", "Debugging trace", TOOL_NAME), false);
//                }
//            } else if (a.equals("-S") || a.equalsIgnoreCase("--stats-file")) {
//                if (it.hasNext()) {
//                    STATS_FILE = it.next();
//                }
//            } else if (a.equals("-N") || a.equalsIgnoreCase("--name-prefix")) {
//                if (it.hasNext()) {
//                    NAME_PREFIX = it.next();
//                }
//            } else if (a.equals("-F") || a.equalsIgnoreCase("--fasta-out")) {
//                OUTPUT_FASTA = true;
//            } else if (a.equals("-O") || a.equalsIgnoreCase("--out-redirect")) {
//                if (it.hasNext()) {
//                    OUT_REDIRECT = it.next();
//                }
//            } else if (a.equals("-E") || a.equalsIgnoreCase("--err-redirect")) {
//                if (it.hasNext()) {
//                    ERR_REDIRECT = it.next();
//                }
//            } else if (a.equalsIgnoreCase("--ideas")) {
//                RUN_SOME_WILD_AND_WONDERFUL_STUFF = true;
////            } else if (a.equals("-M") || a.equalsIgnoreCase("--multipass-compress")) {
////                if (it.hasNext()) {
////                    String optionValue = it.next();
////                    try {
////                        MULTIPASS_COMPRESS = Integer.parseInt(optionValue);
////                    } catch (NumberFormatException e) {
////                        System.err.println("Fatal error! The number iterations set with -m, --multipass-compress must be given as an integer, offending argument: \"" + optionValue + "\"");
////                    }
////                }
////            } else if (a.equals("-H") || a.equalsIgnoreCase("--hash-array-size")) {
////                if (it.hasNext()) {
////                    String optionValue = it.next();
////                    try {
////                        HASH_ARRAY_SIZE = Integer.parseInt(optionValue.replaceAll(",", ""));
////                        if (HASH_ARRAY_SIZE > 4) {
////                            System.err.println("Fatal error! The -H, --hash-array-size must be equal or less than " + (Integer.MAX_VALUE - 3));
////                            System.exit(1);
////                        }
////                    } catch (NumberFormatException e) {
////                        System.err.println("Fatal error! The -H, --hash-array-size must  be given as an integer, offending argument: \"" + optionValue + "\"");
////                    }
////                }
//            } else if (a.equals("-f") || a.equalsIgnoreCase("--min-kmer-frequency")) {
//                if (it.hasNext()) {
//                    String optionValue = it.next();
//                    try {
//                        MIN_KMER_FREQUENCY = Integer.parseInt(optionValue.replaceAll(",", ""));
//                        if (MIN_KMER_FREQUENCY < 1) {
//                            System.err.println("Fatal error! The -f, --min-kmer-frequency must be equal or less than 1");
//                            System.exit(1);
//                        }
//                    } catch (NumberFormatException e) {
//                        System.err.println("Fatal error! The -f, --min-kmer-frequency must  be given as an integer, offending argument: \"" + optionValue + "\"");
//                    }
//                }
//            } else if (a.startsWith("-")) {
//                System.out.println("Unrecognized argument " + a + ", try -h or --help. Terminating...");
//                System.exit(1);
//            } else {
//                inputFileNamesList.add(a);
//            }
//        }
////        if (inputFileNamesList.isEmpty() && MULTIPASS_COMPRESS != null) {
////            System.out.println("Use of -M, --multi-pass-compress requires that an input file is specified, it will not work correctly with stdin");
////            System.exit(1);
////        }
//    }
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
