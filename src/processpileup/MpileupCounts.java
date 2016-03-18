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
import argparser.PositionalOpt;
import gbssplit.SampleBuffer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
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
import java.util.zip.GZIPInputStream;
import kmerextender.KmerExtender;
import kmerextender.PairMerMapPopulatorConsumer;
import shared.Info;
import shared.InputReaderProducer;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class MpileupCounts {

    private final String TOOL_NAME;
    boolean IGNORE_SPLICING = true;
//    private String INPUT_FILE;

//    private ArrayList<String> SAMPLE_NAMES;
    private final int HELP_WIDTH = 200;

    public MpileupCounts(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        parallelMpileupProcessing(optSet);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Process multi-sample mpileup input. Base call lines printed to stdout prefixed with ^CALLS\\t, "
            + "base count lines printed to stdout prefixed by ^COUNTS\\t");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('n', "sample-names", "space separated sample names, order must correspond to mpileup input").setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('p', "pileup-file", "Input (m)pileup file, alternatively use stdin", 1));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Base calling settings]");
        optSet.addOpt(new Opt('a', "min-coverage-per-allele", "Minimum coverage required for an allele to be considered", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt('c', "min-coverage-per-locus", "Minimum coverage required for a locus to be considered", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt('C', "max-coverage-per-locus", "Maximum coverage allowed for a locus to be considered", 1).setMinValue(1).setDefaultValue(1000));
        optSet.addOpt(new Opt('s', "min-samples-within-coverage", "Minimum samples within coverage thresholds required to produce ouput", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt('A', "max-percent-error-allele", "Percentage of coverage up to which an allele is regarded to be an error", 1).setMinValue(0.0).setDefaultValue(1.0));
        optSet.addOpt(new Opt('L', "max-percent-error-locus", "Percentage coverage of alternative alleles up to which a locus is reported", 1).setMinValue(0.0).setDefaultValue(1.0));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Runtime settings]");
        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 1024, 128, 8192));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
            64, 1, 256));
//        String headerNote = "Can be useful for external parallization (print header once)";
//        optSet.addOpt(new Opt('H', "header-only", "Print header and exit").addFootnote(1, TOOL_NAME));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[A little bit of help]");
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
//        System.err.print("Ref\tPos");
//        for (int i = 0; i < sampleIds.size(); i++) {
//            System.err.print("\tA,C,G,T");            
//        }
//        System.err.println("");
        return header;
    }

    private void parallelMpileupProcessing(OptSet optSet) {
        boolean exit = false;
        if (optSet.getOpt("I").isUsed()) {
            System.out.println(Info.getIupacCodesTable());
            exit = true;
        }
        if (optSet.getOpt("D").isUsed()) {
            System.out.println(getOutputCodes());
            exit = true;
        }
        if (exit) {
            System.exit(0);
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
        int MAX_THREADS = (int) optSet.getOpt("t").getValueOrDefault();

//        mpileupCounts(PILEUP_FILE, minCoverageThreshold, maxCoverageThreshold, minSamples, maxPercAlternative);
        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = MAX_THREADS;
            Reporter.report("[INFO]", "Allocated " + threads + " thread(s) to mpileup processing", TOOL_NAME);
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateExecutor = new ThreadPoolExecutor(threads + 1, threads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);

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

            //SPAWN CONSUMER THREADS 
            for (int i = 0; i < threads; i++) {
//                MpileupConsumer consumer = new MpileupConsumer(inputQueue, minCoveragePerLocus, maxCoveragePerLocus, minSamples, maxPercErrAllele, TOOL_NAME);
                MpileupConsumer consumer = new MpileupConsumer(inputQueue, optSet, TOOL_NAME);
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

    private CharSequence getOutputCodes() {
        StringBuilder codes = new StringBuilder();
        codes.append("|---------+--------------------------------------------------|\n");
        codes.append("|  Symbol | Description                                      |\n");
        codes.append("|---------+--------------------------------------------------|\n");
        codes.append("|    ?    + insufficient or conflicting information to call  |\n");
        codes.append("|    d    + deletion in read / insertion in reference        |\n");
        codes.append("|    i    + insertion in read / deletion in the reference    |\n");
        codes.append("|         + whitespace indicates no aligned bases            |\n");
        codes.append("|---------+--------------------------------------------------|\n\n");

        codes.append("|------------------------------------------------------------|\n");
        codes.append("|  COUNTS array represents numbers of identified bases       |\n");
        codes.append("|  [0] ignored and not printed (for now)                     |\n");
        codes.append("|------------------------------------------------------------|\n");
        codes.append("|  [A,C,T,G,d,i]                                             |\n");
        codes.append("|------------------------------------------------------------|\n\n");

        return codes;
    }
}