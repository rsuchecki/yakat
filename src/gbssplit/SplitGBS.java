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
import java.io.File;
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
    private final int MIN_LENGTH_READ;
    private final int MIN_LENGTH_PAIR;
    private final String R1_SUFFIX;
    private final String R2_SUFFIX;
    private final String SE_SUFFIX;
    private final String BLANK_SAMPLE_NAME; //can be repeated so name will get extended with remaining key file columns
    
    private final int HELP_WIDTH = 175;

    public SplitGBS(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        KEY_FILE_NAME = (String) optSet.getOpt("-K").getValueIfSingle();
        BLANK_SAMPLE_NAME = (String) optSet.getOpt("B").getValueOrDefault();
        OUT_DIR = (String) optSet.getOpt("o").getValueOrDefault();        
        new File(OUT_DIR).mkdirs(); //make sure outdir exists
        if (optSet.getOpt("keep-barcodes").getOptFlag()) {
            TRIM_BARCODE = false;
        }
        MIN_LENGTH_READ = (int) optSet.getOpt("r").getValueOrDefault();
        MIN_LENGTH_PAIR = (int) optSet.getOpt("p").getValueOrDefault();
        SPLITTER_THREADS = (int) optSet.getOpt("t").getValueOrDefault();
        R1_SUFFIX = (String) optSet.getOpt("-x").getValueOrDefault();
        R2_SUFFIX = (String) optSet.getOpt("-X").getValueOrDefault();
        SE_SUFFIX = (String) optSet.getOpt("-S").getValueOrDefault();
        
//        for(Opt o: optSet.getOptsList()) {
//            Reporter.report("[INFO]", o.getOptLabelString()+" "+o.getValueOrDefault(), toolName);
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
        //INPUT
        optSet.setListingGroupLabel("[Input options]");
        optSet.addOpt(new Opt('K', "key-file", "Key file name ", 1).setRequired(true));
        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output options]");
        optSet.addOpt(new Opt('o', "out-dir", "Output directory", 1).setDefaultValue("out_split"));
        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
        optSet.addOpt(new Opt('r', "min-length-per-read", "Only output a read if length is no less than <arg> bp", 1, 1, null, 1, 1));
        optSet.addOpt(new Opt('p', "min-length-per-pair", "Only output a read pair if combined length is no less than <arg> bp", 2, 2, null, 1, 1));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime options]");
        optSet.addOpt(new Opt('t', "splitter-threads", "Number of splitter threads. No point setting too high, "
                + "i/o is the likely bottleneck and a writing thread will be spawned per each sample", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        int footId = 1;
        String footText = "Carefully increase to sacrifice memory for speed. Decrease if incountering Out of memory errors.";
        optSet.addOpt(new Opt('u', "buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "handled by a thread (one per sample)", 1024, 64, 8092).addFootnote(footId, footText));
        optSet.addOpt(new Opt('q', "queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up. ", 
                2, 1, 32).addFootnote(footId, footText));
        
        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void splitFiles(ArrayList<String> inputFilenamesList) {
        KeyMap keyMap = new KeyMap(KEY_FILE_NAME, TOOL_NAME, BLANK_SAMPLE_NAME);

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

        BlockingQueue inputQueue = new ArrayBlockingQueue(2);
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
//            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(OUT_DIR + "/" + sample, queue, TOOL_NAME, SPLITTER_THREADS)));
            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(queue, TOOL_NAME, OUT_DIR, sample, SPLITTER_THREADS, R1_SUFFIX, R2_SUFFIX, SE_SUFFIX)));

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
