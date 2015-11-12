/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kmersetmerge;

import agrparser.ArgParser;
import agrparser.Opt;
import agrparser.OptSet;
import agrparser.PositionalOpt;
import agrparser.UsageAndHelp;
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
import shared.InputReaderProducer;
import shared.Reporter;
import shared.WriterConsumer;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerSetMerge {

    private String OUT_LABEL;
    private Integer MIN_FREQUENCY;
    private Integer MAX_FREQUENCY;
    private final String TOOL_NAME;

    public KmerSetMerge(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, 180);
        TOOL_NAME = callerName+" "+toolName;
        try {
            OUT_LABEL = (String) optSet.getOpt("-L").getValueIfSingle();
        } catch (NullPointerException e) {
        }
        try {
            MIN_FREQUENCY = Integer.valueOf((String) optSet.getOpt("-m").getValueIfSingle());
        } catch (NullPointerException | NumberFormatException e) {
        }
        try {
            MAX_FREQUENCY = Integer.valueOf((String) optSet.getOpt("-M").getValueIfSingle());
        } catch (NullPointerException | NumberFormatException e) {
        }
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        mergeSets(inputFilenamesList);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        optSet.addOpt(new Opt('m', "min-frequency", "Output k-mers whose frequency is greater or equal <arg>", 1, 1, Long.MAX_VALUE));
        optSet.addOpt(new Opt('M', "max-frequency", "Output k-mers whose frequency is less then or equal <arg>", null, 1, Long.MAX_VALUE));
        optSet.addOpt(new Opt('L', "re-label", "Set the output k-mers labels to <arg>", 1));
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, Short.MAX_VALUE));
        return optSet;
    }

    private void mergeSets(ArrayList<String> inputFilenamesList) {
        Integer k = -1;
        //READ INPUT AND POPULATE PairMers MAP
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);
        int threads = 3;
        BlockingQueue inputQueue = new ArrayBlockingQueue(65536);
        BlockingQueue outputQueue = new ArrayBlockingQueue(65536);
//            boolean stranded = false;
        ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
        final ExecutorService executorService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFilenamesList, k, TOOL_NAME);
        Future<?> future = executorService.submit(inputReaderProducer);
        futures.add(future);

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
        //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
        if (!inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
            Reporter.report("[FATAL]", "Only k-mer sets accepted as input. Guessed format: " + inputReaderProducer.getGuessedInputFormat(), getClass().getSimpleName());
        } else {
            //Start KmergerConsumerProducer and OutputWriterConsumer threads
            KmergerConsumer kmerger = new KmergerConsumer(inputQueue, outputQueue, MIN_FREQUENCY, MAX_FREQUENCY, OUT_LABEL);
            futures.add(executorService.submit(kmerger));
            WriterConsumer writer = new WriterConsumer(outputQueue, TOOL_NAME);
            futures.add(executorService.submit(writer));
        }
        executorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "PairMerSet populator execution exception!", getClass().getSimpleName());
            ex.getMessage();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "PairMerSet populator timeout exception!", getClass().getSimpleName());
        }
    }

}
