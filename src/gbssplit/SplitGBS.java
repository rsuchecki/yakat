/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbssplit;

import kmerger.*;
import agrparser.ArgParser;
import agrparser.Opt;
import agrparser.OptSet;
import agrparser.PositionalOpt;
import agrparser.UsageAndHelp;
import java.util.ArrayList;
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
import shared.FileWriterConsumer;
import shared.InputReaderProducer;
import shared.Reporter;
import shared.WriterConsumer;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SplitGBS {

    private String KEY_FILE_NAME;
    private String OUT_DIR;
    private final String TOOL_NAME;
    private boolean TRIM_BARCODE = true;

    public SplitGBS(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, 180);
        TOOL_NAME = callerName + " " + toolName;
        try {
            KEY_FILE_NAME = (String) optSet.getOpt("-K").getValueIfSingle();
        } catch (NullPointerException e) {
        }
        OUT_DIR = (String) optSet.getOpt("o").getDefaultValue();
        try {
            String outDir = (String) optSet.getOpt("o").getValueIfSingle();
            if(outDir != null) {
                OUT_DIR = outDir;
            }
        } catch (NullPointerException e) {
        }
//        try {
//            MIN_FREQUENCY = Integer.valueOf((String) optSet.getOpt("-m").getValueIfSingle());
//        } catch (NullPointerException | NumberFormatException e) {
//        }
//        try {
//            MAX_FREQUENCY = Integer.valueOf((String) optSet.getOpt("-M").getValueIfSingle());
//        } catch (NullPointerException | NumberFormatException e) {
//        }
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        splitFiles(inputFilenamesList);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        optSet.addOpt(new Opt('K', "key-file", "Key file name ", 1));
        Opt opt = new Opt('o', "out-dir", "Output directory", 1);
        opt.setDefaultValue("out_split");
        optSet.addOpt(opt);
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, Short.MAX_VALUE));
        return optSet;
    }

    private void splitFiles(ArrayList<String> inputFilenamesList) {
        KeyMap keyMap = new KeyMap(KEY_FILE_NAME, TOOL_NAME);

//        ArrayList<String> r1 = new ArrayList<>();
//        ArrayList<String> r2 = new ArrayList<>();
//        ArrayList<String> se = new ArrayList<>();
//        
//        for (String name : inputFilenamesList) {
//            String tmpName = name.toLowerCase().replace(".gz", "");
//            if(tmpName.endsWith("r1.fastq")) {
//                r1.add(name);
//            } else if(tmpName.endsWith("r2.fastq")) {
//                r2.add(name);
//            } else {
//                se.add(name);
//            }
//        }
//        
//        System.err.println("R1 files:");
//        for (String r : r1) {
//            System.err.println("\t"+r);
//        }
//        System.err.println("R2 files:");
//        for (String r : r2) {
//            System.err.println("\t"+r);
//        }
//        System.err.println("SE files:");
//        for (String r : se) {
//            System.err.println("\t"+r);            
//        }
//        
        Integer k = -1;
        //READ INPUT AND POPULATE PairMers MAP
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);
        int threads = keyMap.getSamplesTotal() + 2;
        BlockingQueue inputQueue = new ArrayBlockingQueue(255);
//            boolean stranded = false;
        ArrayList<Future<?>> futures = new ArrayList<>(threads + 5);
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
        int splitterThreads =1; 
        for (int i = 0; i < splitterThreads; i++) {
            futures.add(executorService.submit(new SplitterConsumerProducer(inputQueue, keyMap, TOOL_NAME, TRIM_BARCODE, splitterThreads)));
        }

        for (Map.Entry<String, BlockingQueue<ArrayList<String>>> entrySet : keyMap.getSampleToQueueMap().entrySet()) {
            String sample = entrySet.getKey();
            BlockingQueue<ArrayList<String>> queue = entrySet.getValue();
            futures.add(executorService.submit(new FileWriterConsumer(OUT_DIR + "/" + sample, queue, TOOL_NAME, splitterThreads)));
//            Reporter.report("[INFO]", "Created writer thread for "+sample, TOOL_NAME);

        }

        executorService.shutdown();
        try {
            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause(), getClass().getSimpleName());
            ex.printStackTrace();

        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
    }

}
