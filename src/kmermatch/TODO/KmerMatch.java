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
package kmermatch.TODO;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptGroup;
import argparser.OptSet;
import argparser.PositionalOpt;
import shared.Reporter;
import shared.InputReaderProducer;
import kmerextender.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerMatch {

    private KmerMap kmersMap = new KmerMap();
    private final ArrayList<String> inputFileNamesList = new ArrayList<>();
    private Integer KMER_LENGTH;
    private Integer MIN_KMER_FREQUENCY;

    private int MAX_THREADS = Runtime.getRuntime().availableProcessors();
    private String DEBUG_FILE;
    private String STATS_FILE;
    private boolean OUTPUT_FASTA = false;
    private String NAME_PREFIX = "";
    private String OUT_REDIRECT;
    private String ERR_REDIRECT;
    private final String TOOL_NAME;

//    private Integer HASH_ARRAY_SIZE = Integer.MAX_VALUE - 3;
//    private Integer MULTIPASS_COMPRESS = null;
    private boolean RUN_SOME_WILD_AND_WONDERFUL_STUFF = false;
    private final int HELP_WIDTH = 180;

    private enum InputType {

        KMERS, FASTA, FASTQ;
    }

    public KmerMatch(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;
        
        
        
        
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        
        System.exit(0);
        
//        processArgs(args);
        
        
        
        
        if (RUN_SOME_WILD_AND_WONDERFUL_STUFF) {

        } else {
            if (OUT_REDIRECT != null) {
                try {
                    File file = new File(OUT_REDIRECT);
                    PrintStream printStream;
                    printStream = new PrintStream(new FileOutputStream(file));
                    System.setOut(printStream);
                } catch (FileNotFoundException ex) {
                    Reporter.report("[ERROR]", "Failed redirecting stdout to " + OUT_REDIRECT, getClass().getSimpleName());
                }
            }
            if (ERR_REDIRECT != null) {
                try {
                    File file = new File(ERR_REDIRECT);
                    PrintStream printStream;
                    printStream = new PrintStream(new FileOutputStream(file));
                    System.setErr(printStream);
                } catch (FileNotFoundException ex) {
                    Reporter.report("[ERROR]", "Failed redirecting stderr to " + ERR_REDIRECT, getClass().getSimpleName());
                }
            }
            runKmerMatcher();
        }
    }
    
    
    
    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Target input files]");
        optSet.addOpt(new Opt('I', "in-single", "Target SE FASTA/FASTQ file(s)", 1).setMaxValueArgs(Short.MAX_VALUE));
        optSet.addOpt(new Opt('1', "in-paired-r1", "Target R1 FASTA/FASTQ file(s), use in conjunction with -2 <args>", 1).setMaxValueArgs(Short.MAX_VALUE));
        optSet.addOpt(new Opt('2', "in-paired-r2", "Target R2 FASTA/FASTQ file(s), use in conjunction with -1 <args>", 1).setMaxValueArgs(Short.MAX_VALUE));
        
        //REFERENCE INPUT TYPE 1 
        int refGroupI = optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel(refGroupI, "[Reference input - option I]");
        optSet.addOpt(new Opt('K', "k-mers-set", "A file contining a pre-computed set of reference k-mers, one k-mer per line", 1).setMaxValueArgs(Short.MAX_VALUE));
        
        //REFERENCE INPUT TYPE 2
        int refGroupII = optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel(refGroupII, "[Reference input - option II]");
        optSet.addOpt(new Opt('f', "fast-to-k-mers", "Input FASTA/FASTQ file(s) to be chopped-up into reference set of k-mers", 1).setMaxValueArgs(Short.MAX_VALUE));
        optSet.addOpt(new Opt('k', "k-mer-size", "Specify the size of k if you don't input a set of k-mers",1).setMinValue(4).setMaxValue(1024));
        optSet.setMutuallyExclusiveGroups(refGroupI, refGroupII);
        optSet.setGroupRequired(refGroupI);
        optSet.setGroupRequired(refGroupII);
//        OptGroup optGroup = new OptGroup();
//        optGroup.addOpt(optFast2kmers);
//        optGroup.addOpt(optKmerSize);
//        optSet.addOptGroup(optGroup);
        
//        optSet.addOpt(new Opt('A', "expected-adapter", "Expected adapter sequence (a fragment will do)", 1).setDefaultValue("AGATCGGAA"));
//        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        
        //
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[?????      settings]");
//        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
//        optSet.addOpt(new Opt('a', "keep-adapters", "Do not trim adapters found next to PstI and MspI sites "));
//        optSet.addOpt(new Opt('p', "keep-non-PstI-starting", "Keep reads (or pairs) which do not start with 'barcodeTGCAG'"));
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Matching settings]");
        optSet.addOpt(new Opt('S', "stranded-matching", "Do not reverse-complement k-mers (by default canonical representation of a k-mer is stored and matched)"));
        optSet.addOpt(new Opt('M', "min-matches", "Set minimum number of reference k-mers matching each targeted sequence",1,1,Integer.MAX_VALUE));
        optSet.addOpt(new Opt('V', "invert-matching", "Invert matching: output targets that contain fewer than [-M] k-mers matching the reference set"));
        optSet.addOpt(new Opt('B', "match-both-mates", "Relevant for PE input only. Demand each mate to have [-M] matching kmers with the reference set, by default, both mates are caught if at least one has [-M] kmer(s) matching the reference"));
        optSet.addOpt(new Opt('A', "report-all", "Report all input sequences, with the number of matching baiting k-mers appended to the identifier, other matchin settings will be ignored"));

//        int footId = 1;
//        String footText1 = "Note that certain combinations of min-length-* settings can lead to both mates of a pair ending up in SE/orphans output file.";
//        optSet.addOpt(new Opt('r', "min-length-single-read", "Only output a read if length is no less than <arg> bp", 1, 1, null, 1, 1).addFootnote(footId, footText1));
//        optSet.addOpt(new Opt('e', "min-length-pair-each", "Only output a read pair if length of each is no less than <arg> bp, otherwise process as single", 1, 1, null, 1, 1).addFootnote(footId, footText1));
//        optSet.addOpt(new Opt('s', "min-length-pair-sum", "Only output a read pair if combined length is no less than <arg> bp, otherwise process as single", 2, 2, null, 1, 1).addFootnote(footId, footText1));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
            + "passed to in-queue", 1024, 128, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
            2, 1, 256));
        optSet.addOpt(new Opt('T', "querying-threads", "Number of threads to query the generated reference set. No point setting too high for small reference sets, "
            + "i/o is the likely bottleneck.", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
        optSet.addOpt(new Opt('d', "out-dir", "Output directory", 1).setDefaultValue("matched"));
        optSet.addOpt(new Opt('u', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
        optSet.addOpt(new Opt('x', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
        optSet.addOpt(new Opt('s', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
        int footId = 1;
        String footText2 = "Consider increasing to sacrifice memory for speed. Decrease if encountering 'out of memory' errors.";
        optSet.addOpt(new Opt('r', "out-buffer-size", "Number of FASTQ records (reads or pairs) "
            + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText2));
        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
            2, 1, 32).addFootnote(footId, footText2));

        //POSITIONAL
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }
    
    
    
    

    private void runKmerMatcher() {
        Reporter.report("[INFO]", "Initialized, will use " + MAX_THREADS + " thread(s) to populate map ", getClass().getSimpleName());
        readKmersAndPopulateKmersMap();

        Reporter.report("[INFO]", "Finished populating map, k=" + KMER_LENGTH + ", n=" + NumberFormat.getIntegerInstance().format(kmersMap.getPairMersSkipListMap().size()), getClass().getSimpleName());
//        PairMersExtender pairMersExtender = new PairMersExtender(DEBUG_FILE, STATS_FILE);
//        pairMersExtender.matchAndExtendKmers(KMER_LENGTH, kmersMap, OUTPUT_FASTA, NAME_PREFIX, MAX_THREADS);
//        Reporter.report("[INFO]", "Finished extending k-mers");
    }

    /**
     * Producer-consumer multi-threading to read input (k-mers, FASTA, etc)
     * generate Kmer objects and populate a Map
     */
    private void readKmersAndPopulateKmersMap() {
        try {
            Integer k;
            //READ INPUT AND POPULATE PairMers MAP
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);
            int threads = MAX_THREADS;
            BlockingQueue inputQueue = new ArrayBlockingQueue(65536);
//            boolean stranded = false;
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateDbExecutor = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            //SPAWN INPUT READING THREAD
            InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFileNamesList, KMER_LENGTH, kmersMap, TOOL_NAME);
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
                        Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", getClass().getSimpleName());
                    }
                    Thread.sleep(100); //wait for 1/10 of a second
                } catch (InterruptedException ex) {
                }
            }
            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
                splitInputSequenceIntoKmers = false;
                KMER_LENGTH = inputReaderProducer.getKmerLength();
            }

            //SPAWN THREADS TO POPULATE CLIPMERS MAP
            for (int i = 0; i < threads; i++) {
                KmerMapPopulatorConsumer consumer = new KmerMapPopulatorConsumer(inputQueue, kmersMap, splitInputSequenceIntoKmers, KMER_LENGTH, MIN_KMER_FREQUENCY);
                futures.add(readAndPopulateDbExecutor.submit(consumer));
            }
            readAndPopulateDbExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                readAndPopulateDbExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!", getClass().getSimpleName());
            } catch (ExecutionException ex) {
                Reporter.report("[ERROR]", "PairMerSet populator execution exception!", getClass().getSimpleName());
                ex.getMessage();

            } catch (TimeoutException ex) {
                Logger.getLogger(KmerMatch.class.getName()).log(Level.SEVERE, null, ex);
                Reporter.report("[ERROR]", "PairMerSet populator timeout exception!", getClass().getSimpleName());
            }

            if (kmersMap.isOutOfMemory()) {
                Reporter.report("[ERROR]", "Terminating, out of memory while populating the Map, k=" + KMER_LENGTH + ", n=" + NumberFormat.getIntegerInstance().format(kmersMap.getPairMersSkipListMap().size()), getClass().getSimpleName());
                System.exit(1);
            }

        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error!", getClass().getSimpleName());
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            System.exit(1);
        }

    }


    /**
     * Supposedly POSIX-compliant CLI args processor
     *
     * @param args
     */
    private void processArgs(String[] args) {
        ArrayList<String> allArgs = new ArrayList<>();
        //SPLIT-UP POSIX STYLE (-abcd) IF ANY
        for (String arg : args) {
            if (arg.startsWith("--")) {
                allArgs.add(arg);
            } else if (arg.startsWith("-")) {
                if (arg.length() > 2) {
                    char[] chars = arg.toCharArray();
                    StringBuilder sb = new StringBuilder();
                    for (int i = 1; i < chars.length; i++) {
                        char c = chars[i];
//                        System.out.println(c + " ASCII value " + (int) c);
                        if ((int) c >= 48 && (int) c <= 57) {
                            sb.append(c);
                        } else {
                            allArgs.add("-" + c);
                        }
                    }
                    String toString = sb.toString();
                    if (!toString.isEmpty()) {
                        allArgs.add(toString);
                    }

                } else {
                    allArgs.add(arg);
                }
            } else {
                allArgs.add(arg);
            }
        }
        //Now process args
        Iterator<String> it = allArgs.iterator();

        while (it.hasNext()) {
            String a = it.next();
            //
            if (a.equals("-h") || a.equals("--help")) {
                printHelp();
                System.exit(0);
            } else if (a.equals("-k") || a.equalsIgnoreCase("--kmer-length")) {
                if (it.hasNext()) {
                    String optionValue = it.next();
                    try {
                        KMER_LENGTH = Integer.parseInt(optionValue);
                        if (KMER_LENGTH < 4) {
                            System.err.println("Fatal error! \"k\" must be set to at least 4, currently set to: " + KMER_LENGTH);
                            System.exit(1);
                        }
                    } catch (NumberFormatException e) {
                        System.err.println("Fatal error! The k-mer lenght must be given as an integer, offending argument: \"" + optionValue + "\"");
                    }
                }
            } else if (a.equals("-t") || a.equalsIgnoreCase("--threads")) {
                if (it.hasNext()) {
                    String optionValue = it.next();
                    try {
                        MAX_THREADS = Integer.parseInt(optionValue);
                        if (MAX_THREADS < 1) {
                            System.err.println("Fatal error! \"threads\" must be set to a positive integer, currently set to: " + MAX_THREADS);
                            System.exit(1);
                        } else if (MAX_THREADS > Runtime.getRuntime().availableProcessors()) {
                            MAX_THREADS = Runtime.getRuntime().availableProcessors();
                            System.err.println("Warning, number of threads reset to the maximum available (" + MAX_THREADS + ") CPU cores, offending argument: \"" + optionValue + "\"");
                        }
                    } catch (NumberFormatException e) {
                        System.err.println("Fatal error! Number of threads must be given as an integer, offending argument: \"" + optionValue + "\"");
                    }
                }
            } else if (a.equals("-D") || a.equalsIgnoreCase("--debug-file")) {
                if (it.hasNext()) {
                    DEBUG_FILE = it.next();
                    Reporter.writeToFile(DEBUG_FILE, Reporter.formatReport("[DEBUG]", "Debugging trace", getClass().getSimpleName()), false);
                }
            } else if (a.equals("-S") || a.equalsIgnoreCase("--stats-file")) {
                if (it.hasNext()) {
                    STATS_FILE = it.next();
                }
            } else if (a.equals("-N") || a.equalsIgnoreCase("--name-prefix")) {
                if (it.hasNext()) {
                    NAME_PREFIX = it.next();
                }
            } else if (a.equals("-F") || a.equalsIgnoreCase("--fasta-out")) {
                OUTPUT_FASTA = true;
            } else if (a.equals("-O") || a.equalsIgnoreCase("--out-redirect")) {
                if (it.hasNext()) {
                    OUT_REDIRECT = it.next();
                }
            } else if (a.equals("-E") || a.equalsIgnoreCase("--err-redirect")) {
                if (it.hasNext()) {
                    ERR_REDIRECT = it.next();
                }
            } else if (a.equalsIgnoreCase("--ideas")) {
                RUN_SOME_WILD_AND_WONDERFUL_STUFF = true;
//            } else if (a.equals("-M") || a.equalsIgnoreCase("--multipass-compress")) {
//                if (it.hasNext()) {
//                    String optionValue = it.next();
//                    try {
//                        MULTIPASS_COMPRESS = Integer.parseInt(optionValue);
//                    } catch (NumberFormatException e) {
//                        System.err.println("Fatal error! The number iterations set with -m, --multipass-compress must be given as an integer, offending argument: \"" + optionValue + "\"");
//                    }
//                }
//            } else if (a.equals("-H") || a.equalsIgnoreCase("--hash-array-size")) {
//                if (it.hasNext()) {
//                    String optionValue = it.next();
//                    try {
//                        HASH_ARRAY_SIZE = Integer.parseInt(optionValue.replaceAll(",", ""));
//                        if (HASH_ARRAY_SIZE > 4) {
//                            System.err.println("Fatal error! The -H, --hash-array-size must be equal or less than " + (Integer.MAX_VALUE - 3));
//                            System.exit(1);
//                        }
//                    } catch (NumberFormatException e) {
//                        System.err.println("Fatal error! The -H, --hash-array-size must  be given as an integer, offending argument: \"" + optionValue + "\"");
//                    }
//                }
            } else if (a.equals("-f") || a.equalsIgnoreCase("--min-kmer-frequency")) {
                if (it.hasNext()) {
                    String optionValue = it.next();
                    try {
                        MIN_KMER_FREQUENCY = Integer.parseInt(optionValue.replaceAll(",", ""));
                        if (MIN_KMER_FREQUENCY < 1) {
                            System.err.println("Fatal error! The -f, --min-kmer-frequency must be equal or less than 1");
                            System.exit(1);
                        }
                    } catch (NumberFormatException e) {
                        System.err.println("Fatal error! The -f, --min-kmer-frequency must  be given as an integer, offending argument: \"" + optionValue + "\"");
                    }
                }
            } else if (a.startsWith("-")) {
                System.out.println("Unrecognized argument " + a + ", try -h or --help. Terminating...");
                System.exit(1);
            } else {
                inputFileNamesList.add(a);
            }
        }
//        if (inputFileNamesList.isEmpty() && MULTIPASS_COMPRESS != null) {
//            System.out.println("Use of -M, --multi-pass-compress requires that an input file is specified, it will not work correctly with stdin");
//            System.exit(1);
//        }
    }

    private void printHelp() {

        System.out.println();
        System.out.println("Given a set of k-mers, performs their unmabiguous extension.");
        System.out.println("Input may be passed through stdin. Output is printed to stdout. ");
        System.out.println("Progress information, warnings errors etc. are printed to stderr, usage:");
        System.out.println("java -jar " + this.getClass().getSimpleName() + ".jar [options] [INPUT_FILENAME(s)] ");
        System.out.println(" -h, --help                     : Print this help screen and exit");
        System.out.println(" -k, --kmer-length <int>        : Required only if input other than a list of k-mers");
        System.out.println(" -t, --threads <int>            : Number of threads, defaults to " + MAX_THREADS);
        System.out.println("                                : ");
        System.out.println(" -F, --fasta-out                : Output each k-mer as a separate FASTA record instead of just listing extended nucleotide sequences");
        System.out.println(" -N, --name-prefix <string>     : If outputing FASTA, prefix each record identifier with this <string>");
        System.out.println("                                : ");
        System.out.println(" -O, --out-redirect <string>    : Redirect stdout to this file");
        System.out.println(" -E, --err-redirect <string>    : Redirect stderr to this file");
        System.out.println(" -S, --stats-file <string>      : Write extension stats to this file");
        System.out.println(" -D, --debug-file <string>      : Write unkosher extensions details to this file");
        System.out.println("                                : ");
//        System.out.println(" Multi-pass hashing options, incompatible with stdin:");
//        System.out.println(" -M, --multipass-compress <int> : Reduce memory requirements with <int> hashing iterations");
//        System.out.println(" -H, --hash-array-size    <int> : At each iteration, array of <int> size will be added, defualts to " + NumberFormat.getNumberInstance().format(HASH_ARRAY_SIZE) + "*");
//        System.out.println("                                : ");
//
//        String s = "* - Note that each hash array will consume 4 Bytes per element, so the largest possible (on most VMs) "
//                + "2,147,483,644 element array will consume over 8GB of the allocated memory. ";
//        System.out.println(Reporter.wrapString(s, 145));
        String s = "Currently k-mer frequency is not taken into consideration, so use of a dedicated k-mer counting program, "
                + "such as KMC or Jellyfish is recommended. It is best to exclude low frequency k-mers before passing "
                + "the list of k-mers to KmerExtender. For smaller jobs FASTA or FASTQ input may suffice.";
        System.out.println(Reporter.wrapString(s, 145));
        System.out.println();
    }

}
