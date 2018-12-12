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
package processpileup;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import java.io.PrintStream;
import java.util.ArrayList;
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
import kextender.KmerExtender;
import shared.Info;
import shared.InputReaderProducer;
import shared.Reporter;
import shared.StdRedirect;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class ProcessPileup {

    private final String TOOL_NAME;
//    boolean IGNORE_SPLICING = true;
//    private String INPUT_FILE;

//    private ArrayList<String> SAMPLE_NAMES;
    private final int HELP_WIDTH = 200;

    public ProcessPileup(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
//        stdRedirect(optSet);
        new StdRedirect(optSet, TOOL_NAME);
        parallelMpileupProcessing(optSet);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Process multi-sample mpileup input. Base calls printed to stdout on lines prefixed with ^CALLS\\t, "
            + "base count printed to stdout on lines prefixed by ^COUNTS\\t");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('n', "sample-names", "space separated sample names, order must correspond to mpileup input").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('p', "pileup-file", "Input (m)pileup file, alternatively use stdin", 1));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Base calling settings]");
        optSet.addOpt(new Opt('a', "min-coverage-per-allele", "Minimum coverage required for an allele to be considered in a locus call", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt('c', "min-coverage-per-locus", "Minimum coverage required for a locus to be considered ", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt('C', "max-coverage-per-locus", "Maximum coverage allowed for a locus to be considered", 1).setMinValue(1).setDefaultValue(1000));
        optSet.addOpt(new Opt('A', "max-percent-error-allele", "Percentage of coverage up to which an allele is regarded to be an error", 1).setMinValue(0.0).setDefaultValue(1.0));
        optSet.addOpt(new Opt('L', "max-percent-error-locus", "Percentage coverage of alternative alleles up to which a locus is reported", 1).setMinValue(0.0).setDefaultValue(1.0));
        optSet.addOpt(new Opt(null, "min-minor-major-ratio", "[TODO] Minimum fraction of a minor allele bases required to call a heterozygous base", 1).setMinValue(0.0).setMaxValue(0.5).setDefaultValue(0.0));
        optSet.addOpt(new Opt(null, "zero-reads-char", "A character denoting zero reads at a postion for a given sample", 1).setDefaultValue('.'));
        optSet.addOpt(new Opt(null, "ambiguous-call-char", "A character indicating uncertain call (e.g. due to low coverage or unclear zygosity at locus)", 1).setDefaultValue('?'));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Record reporting settings]");
        optSet.addOpt(new Opt('s', "min-samples-within-coverage", "Minimum samples within coverage thresholds required to produce ouput", 1).setMinValue(0).setDefaultValue(2));
        
        optSet.addOpt(new Opt('W', "all-within-thresholds", "Print all loci within given thresholds even if no alternative alleles called"));
        optSet.addOpt(new Opt(null, "min-samples-called", "[TODO] Minimum samples for which the base was called", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "max-samples-called", "[TODO] Maximum samples for which the base was called", 1).setMinValue(1));
        optSet.addOpt(new Opt(null, "min-samples-zero-coverage", "[TODO] May be useful for presence-absence analyses", 1).setMinValue(0).setDefaultValue(0));
        optSet.addOpt(new Opt(null, "max-samples-zero-coverage", "[TODO] May be useful for presence-absence analyses", 1).setMinValue(0));
        
        optSet.addOpt(new Opt(null, "min-snps-to-reference", "Report a position if at least <arg> samples have a SNP to reference ",1 ).setMinValue(1));
        optSet.addOpt(new Opt(null, "min-calls-uncertain", "Minimum samples for which the base was not called due to uncertainty", 1).setMinValue(0).setDefaultValue(0));
        optSet.addOpt(new Opt(null, "max-calls-uncertain", "Maximum samples for which the base was not called due to uncertainty", 1).setMinValue(0).setDefaultValue(65535));
        optSet.addOpt(new Opt(null, "min-calls-het", "Report a position if at least <arg> calls are heterozygous",1 ).setMinValue(0).setDefaultValue(0));
        optSet.addOpt(new Opt(null, "max-calls-het", "Report a position if at most <arg> calls are heterozygous",1 ).setMinValue(0).setDefaultValue(65535));
//        optSet.addOpt(new Opt('z', "min-missing-samples", "Minimum samples with zero coverage", 1).setMinValue(0).setDefaultValue(0));
//        optSet.addOpt(new Opt('u', "max-uncalled-samples", "Maximum samples for which the base was not called", 1).setMinValue(0).setDefaultValue(0));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Runtime settings]");
        String threadsOrderNote = "Note that in multi-threaded mode the output lines order need not reflect the input order";
        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(1, threadsOrderNote));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 1024, 128, 32768));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",64, 1, 256));
//        optSet.addOpt(new Opt('u', "out-buffer-size", "Size of buffers put on out-queue ", 1024, 128, 32768));
//        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writing-out",64, 1, 256));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
//        String headerNote = "Can be useful for external parallization (print header once)";
//        optSet.addOpt(new Opt('H', "header-only", "Print header and exit").addFootnote(1, TOOL_NAME));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[A little bit of help]");
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
        optSet.addOpt(new Opt('I', "iupac-codes-table", "Print the table of IUPAC nucleotide codes and exit"));
        optSet.addOpt(new Opt('D', "additional-codes-table", "Print the table of additional codes/symbols used by this program"));
        return optSet;
    }

    private CharSequence getHeader(ArrayList<String> sampleIds) {
        StringBuilder header = new StringBuilder();
        header.append("Ref\tPos\tRefBase");
        for (String sampleId : sampleIds) {
            header.append("\t").append(sampleId);
        }
//        header.append("\t").append("Global");
//        System.err.print("Ref\tPos");
//        for (int i = 0; i < sampleIds.size(); i++) {
//            System.err.print("\tA,C,G,T");            
//        }
//        System.err.println("");
        return header;
    }

//    private void stdRedirect(OptSet optSet) {
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
//    }

    private void parallelMpileupProcessing(OptSet optSet) {

        boolean exit = false;
        if (optSet.getOpt("I").isUsed()) {
            System.out.println(Info.getIupacCodesTable());
            exit = true;
        }
        if (optSet.getOpt("D").isUsed()) {
            System.out.println(getOutputCodes(optSet));
            exit = true;
        }
        if (exit) {
            System.exit(0);
        }
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        ArrayList<String> SAMPLE_NAMES = (ArrayList<String>) optSet.getOpt("n").getValues();
//        int minCoveragePerLocus = (int) optSet.getOpt("c").getValueOrDefault();
//        int maxCoveragePerLocus = (int) optSet.getOpt("C").getValueOrDefault();
//        int minSamples = (int) optSet.getOpt("s").getValueOrDefault();
//        int minCoveragePerAllele = (int) optSet.getOpt("a").getValueOrDefault();
//        
//        double maxPercErrAllele = (double) optSet.getOpt("A").getValueOrDefault();
//        double maxPercErrLocus = (double) optSet.getOpt("L").getValueOrDefault();
        String inputFile = (String) optSet.getOpt("p").getValueOrDefault();
        int READER_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        
        int IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();
//        int OUT_Q_CAPACITY = (int) optSet.getOpt("q").getValueOrDefault();
        int MAX_THREADS = (int) optSet.getOpt("t").getValueOrDefault();

//        mpileupCounts(PILEUP_FILE, minCoverageThreshold, maxCoverageThreshold, minSamples, maxPercAlternative);
        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = MAX_THREADS;
            Reporter.report("[INFO]", "Allocated " + threads + " thread(s) to mpileup processing", TOOL_NAME);
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateExecutor = new ThreadPoolExecutor(threads + 1, threads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
//            BlockingQueue outputQueue = new ArrayBlockingQueue(OUT_Q_CAPACITY);

            //SPAWN INPUT READING THREAD            
            InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFile, null, TOOL_NAME, READER_BUFFER_SIZE);

            Future<?> future = readAndPopulateExecutor.submit(inputReaderProducer);
            futures.add(future);

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

            System.out.println("COUNTS\t" + getHeader(SAMPLE_NAMES));
            System.out.println("CALLS\t" + getHeader(SAMPLE_NAMES));

            PrintStream bufferedOut = new PrintStream(new java.io.BufferedOutputStream(System.out, 65535));
            
            //SPAWN CONSUMER THREADS 
            for (int i = 0; i < threads; i++) {
//                MpileupConsumer consumer = new MpileupConsumer(inputQueue, minCoveragePerLocus, maxCoveragePerLocus, minSamples, maxPercErrAllele, TOOL_NAME);
                MpileupConsumer consumer = new MpileupConsumer(inputQueue, optSet, TOOL_NAME, bufferedOut);
                futures.add(readAndPopulateExecutor.submit(consumer));
            }
            readAndPopulateExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                readAndPopulateExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Reporter.report("[ERROR]", "interrupted exception!", TOOL_NAME);
            } catch (ExecutionException ex) {
                Reporter.report("[ERROR]", "execution exception!", TOOL_NAME);
                ex.printStackTrace();
                System.exit(1);

            } catch (TimeoutException ex) {
                Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
                Reporter.report("[ERROR]", "timeout exception!", TOOL_NAME);
            }

        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error!", TOOL_NAME);
            System.exit(1);
        }
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private CharSequence getOutputCodes(OptSet optSet) {
        StringBuilder codes = new StringBuilder();
        char m = (Character)optSet.getOpt("ambiguous-call-char").getValueOrDefault();
        char a = (Character)optSet.getOpt("zero-reads-char").getValueOrDefault();
        codes.append("|---------+-------------------------------------------------------|\n");
        codes.append("|  Symbol | Description                                           |\n");
        codes.append("|---------+-------------------------------------------------------|\n");
        codes.append("|    "+m+"    + insufficient or conflicting information so not called |\n");
        codes.append("|    E    + deletion in read / insertion in the reference         |\n");
        codes.append("|    I    + insertion in read / deletion in the reference         |\n");
        codes.append("|    "+a+"    + indicates no aligned bases                            |\n");
        codes.append("|---------+-------------------------------------------------------|\n\n");

        codes.append("|------------------------------------------------------------|\n");
        codes.append("|  COUNTS array represents numbers of identified bases       |\n");
        codes.append("|  [0] ignored and not printed (for now)                     |\n");
        codes.append("|------------------------------------------------------------|\n");
        codes.append("|  [A,C,T,G,E,I]                                             |\n");
        codes.append("|------------------------------------------------------------|\n\n");

        return codes;
    }
}
