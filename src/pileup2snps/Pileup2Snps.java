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

package pileup2snps;

import agrparser.*;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import shared.InputReaderProducer;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Pileup2Snps {
    private final int IN_BUFFER_SIZE;
    private final int IN_Q_CAPACITY ;
    
    private final String TOOL_NAME;
    private final int PROCESSING_THREADS = 1;
    private final int HELP_WIDTH = 175;

    public Pileup2Snps(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();
        
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        processPileup(inputFilenamesList);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input options]");
        optSet.addOpt(new Opt('K', "key-file", "Key file name ", 1).setRequired(true));
        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to in-queue", 1024, 128, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up. ", 
                2, 1, 256));
        //OUTPUT & RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime and output options]");
        optSet.addOpt(new Opt('o', "out-dir", "Output directory", 1).setDefaultValue("out_split"));
        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
        optSet.addOpt(new Opt('r', "min-length-per-read", "Only output a read if length is no less than <arg> bp", 1, 1, null, 1, 1));
        optSet.addOpt(new Opt('p', "min-length-per-pair", "Only output a read pair if combined length is no less than <arg> bp", 2, 2, null, 1, 1));
        optSet.addOpt(new Opt('t', "splitter-threads", "Number of splitter threads. No point setting too high, "
                + "i/o is the likely bottleneck and a writing thread will be spawned per each sample", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        int footId = 1;
        String footText = "Consider increasing to sacrifice memory for speed. Decrease if incountering Out of memory errors.";
        optSet.addOpt(new Opt('u', "out-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText));
        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up. ", 
                2, 1, 32).addFootnote(footId, footText));
        
        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void processPileup(ArrayList<String> inputFilenamesList) {
        Integer k = -1;
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);
        ArrayBlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        int ioThreads = 1;
        ArrayList<Future<?>> ioFutures = new ArrayList<>(ioThreads);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(ioThreads, ioThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFilenamesList, k, TOOL_NAME, "records", IN_BUFFER_SIZE);
        ioFutures.add(ioExecutorService.submit(inputReaderProducer));

        long timeStart = System.currentTimeMillis();
        int count = 0;
        while (inputReaderProducer.getGuessedInputFormat() == null) {
            try {
                //IF nothing happens after 5 seconds
                if (System.currentTimeMillis() - timeStart > 2500 && (count++ % 25 == 0)) {
                    Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", getClass().getSimpleName());
                }
                Thread.sleep(100); //wait for 1/10 of a second
            } catch (InterruptedException ex) {
            }
        }
//        //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
//        if (!inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
//            Reporter.report("[FATAL]", "Only k-mer sets accepted as input. Guessed format: " + inputReaderProducer.getGuessedInputFormat(), getClass().getSimpleName());
//        } else {
        //Start KmergerConsumerProducer and OutputWriterConsumer threads
        
        final ExecutorService splitterExecutorService = new ThreadPoolExecutor(PROCESSING_THREADS, PROCESSING_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> splittersFutures = new ArrayList<>(PROCESSING_THREADS);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        for (int i = 0; i < PROCESSING_THREADS; i++) {
            splittersFutures.add(splitterExecutorService.submit(new PileupConsumer(inputQueue, TOOL_NAME, 100)));
        }

//        //WRITER THREADS
//        for (Map.Entry<String, BlockingQueue<ArrayList<String>>> entrySet : keyMap.getSampleToQueueMap().entrySet()) {
//            String sample = entrySet.getKey();
//            BlockingQueue<ArrayList<String>> queue = entrySet.getValue();
////            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(OUT_DIR + "/" + sample, queue, TOOL_NAME, SPLITTER_THREADS)));
//            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(queue, TOOL_NAME, OUT_DIR, sample, SPLITTER_THREADS, R1_SUFFIX, R2_SUFFIX, SE_SUFFIX)));
//
//        }
        splitterExecutorService.shutdown();
        ioExecutorService.shutdown();
        try {

            for (Future<?> f : splittersFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            splitterExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            for (Future<?> f : ioFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            ioExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
    }
}
