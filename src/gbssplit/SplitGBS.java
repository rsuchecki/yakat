/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbssplit;

import agrparser.ArgParser;
import agrparser.Opt;
import agrparser.OptSet;
import agrparser.PositionalOpt;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.atomic.AtomicInteger;
import shared.FileWriterConsumer;
import shared.InputReaderProducer;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SplitGBS {

    private String KEY_FILE_NAME;
    private String OUT_DIR;
    private final String TOOL_NAME;
    private int SPLITTER_THREADS = 1;
    private boolean TRIM_BARCODE = true;
    private Integer MIN_LENGTH_READ;
    private Integer MIN_LENGTH_PAIR;

    public SplitGBS(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, 180);
        TOOL_NAME = callerName + " " + toolName;
        KEY_FILE_NAME = (String) optSet.getOpt("-K").getValueIfSingle();
        OUT_DIR = (String) optSet.getOpt("o").getDefaultValue();
        try {
            String outDir = (String) optSet.getOpt("o").getValueIfSingle();
            if (outDir != null) {
                OUT_DIR = outDir;
            }
        } catch (NullPointerException e) {
        }
        if (optSet.getOpt("keep-barcodes").getOptFlag()) {
            TRIM_BARCODE = false;
        }
        if (optSet.getOpt("r").isUsed()) {
            MIN_LENGTH_READ = Integer.valueOf((String) optSet.getOpt("r").getValueIfSingle());
        }
        if (optSet.getOpt("p").isUsed()) {
            MIN_LENGTH_PAIR = Integer.valueOf((String) optSet.getOpt("p").getValueIfSingle());
        }
        if (optSet.getOpt("t").isUsed()) {
            SPLITTER_THREADS = Integer.valueOf((String) optSet.getOpt("t").getValueIfSingle());
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
        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
        optSet.addOpt(new Opt('r', "min-length-per-read", "Only output a read if length is no less than <arg> bp", null, 1, null, 1, 1));
        optSet.addOpt(new Opt('p', "min-length-per-pair", "Only output a read pair if combined length is no less than <arg> bp", null, 1, null, 1, 1));
        optSet.addOpt(new Opt('t', "splitter-threads", "Number of splitter threads. No point setting too high, "
                + "writing samples out is the likely bottleneck and a separate writing thread will be spawned per sample", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
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

        BlockingQueue inputQueue = new ArrayBlockingQueue(255);
//            boolean stranded = false;
        int ioThreads = keyMap.getSamplesTotal() + 1;
        ArrayList<Future<?>> ioFutures = new ArrayList<>(ioThreads);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(ioThreads, ioThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFilenamesList, k, TOOL_NAME);
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
        final ExecutorService splitterExecutorService = new ThreadPoolExecutor(SPLITTER_THREADS, SPLITTER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> splittersFutures = new ArrayList<>(SPLITTER_THREADS);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        for (int i = 0; i < SPLITTER_THREADS; i++) {
            splittersFutures.add(splitterExecutorService.submit(new SplitterConsumerProducer(inputQueue, keyMap, TOOL_NAME, TRIM_BARCODE, MIN_LENGTH_READ, MIN_LENGTH_PAIR)));
        }
        
        //WRITER THREADS
        for (Map.Entry<String, BlockingQueue<ArrayList<String>>> entrySet : keyMap.getSampleToQueueMap().entrySet()) {
            String sample = entrySet.getKey();
            BlockingQueue<ArrayList<String>> queue = entrySet.getValue();
            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(OUT_DIR + "/" + sample, queue, TOOL_NAME, SPLITTER_THREADS)));

        }
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
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
    }

}
