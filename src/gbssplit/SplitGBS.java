/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gbssplit;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
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
import shared.FileWriterConsumer;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class SplitGBS {

    private final String TOOL_NAME;
    private int SPLITTER_THREADS = 1;

    //IN
    private final int IN_BUFFER_SIZE;
    private final int IN_Q_CAPACITY;

    //KEY
    private final String KEY_FILE_NAME;
    private final String BLANK_SAMPLE_NAME; //can be repeated so name will get extended with remaining key file columns

    //OUT
    private final String OUT_DIR;
    private final String R1_SUFFIX;
    private final String R2_SUFFIX;
    private final String SE_SUFFIX;
    private  String MATCHLESS_OUT = null;
    private final int OUT_BUFFER_SIZE;
    private final int OUT_Q_CAPACITY;
    private boolean OVERWRITE = false;
    private boolean APPEND = false;

    private final int HELP_WIDTH = 180;

    public SplitGBS(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        if (optSet.getOpt("append").isUsed()) {
            APPEND = true;
        }
        if (optSet.getOpt("force").isUsed()) {
            OVERWRITE = true;
        }
        
        if(OVERWRITE && APPEND) {
            Reporter.report("[FATAL]", "Flags --force and --append are mutually exclusive, you can either force-overwrite existing files OR append to them", TOOL_NAME);
            System.exit(1);
        }
        IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();

        BLANK_SAMPLE_NAME = (String) optSet.getOpt("B").getValueOrDefault();
        KEY_FILE_NAME = (String) optSet.getOpt("-K").getValueIfSingle();
        OUT_DIR = (String) optSet.getOpt("o").getValueOrDefault();
        new File(OUT_DIR).mkdirs(); //make sure outdir exists

        R1_SUFFIX = (String) optSet.getOpt("-x").getValueOrDefault();
        R2_SUFFIX = (String) optSet.getOpt("-X").getValueOrDefault();
        SE_SUFFIX = (String) optSet.getOpt("-S").getValueOrDefault();

        SPLITTER_THREADS = (int) optSet.getOpt("t").getValueOrDefault();

        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
        OUT_Q_CAPACITY = (int) optSet.getOpt("q").getValueOrDefault();

        if (optSet.getOpt("M").isUsed()) {
            MATCHLESS_OUT = (String) optSet.getOpt("M").getValueOrDefault();
        }
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }        
        
//        for(Opt o: optSet.getOptsList()) {
//            Reporter.report("[INFO]", o.getOptLabelString()+" "+o.getValueOrDefault(), toolName);
//        }
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        splitFiles(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('K', "key-file", "Key file name ", 1).setRequired(true));
        optSet.addOpt(new Opt(null, "restriction-remnant-I", "First restriction remnant site", 1).setDefaultValue("TGCAG"));
        optSet.addOpt(new Opt(null, "restriction-remnant-II", "Second restriction remnant site", 1).setDefaultValue("CCG"));
        optSet.addOpt(new Opt('A', "expected-adapter", "Expected adapter sequence (a fragment will do)", 1).setDefaultValue("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"));
        optSet.addOpt(new Opt(null, "adapter-prefix-length", "Length of the adapter prefix used to identify 3' read-through", 1).setDefaultValue(9));
        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
            + "passed to in-queue", 1024, 128, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
            2, 1, 256));
        //TRIMMING AND LENGTH
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Trimming and length settings]");
        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
        optSet.addOpt(new Opt('a', "keep-adapters", "Do not trim adapters found next to remnant sites "));
        optSet.addOpt(new Opt('p', "keep-non-remnant-starting", "Keep reads (or pairs) which do not start with barcode followed by first restriction remnant "));

        int footId = 1;
        String footText1 = "Note that certain combinations of min-length-* settings can lead to both mates of a pair ending up in SE/orphans output file.";
        optSet.addOpt(new Opt('r', "min-length-single-read", "Only output a read if length is no less than <arg> bp", 1, 1, null, 1, 1).addFootnote(footId, footText1));
        optSet.addOpt(new Opt('e', "min-length-pair-each", "Only output a read pair if length of each is no less than <arg> bp, otherwise process as single", 1, 1, null, 1, 1).addFootnote(footId, footText1));
        optSet.addOpt(new Opt('s', "min-length-pair-sum", "Only output a read pair if combined length is no less than <arg> bp, otherwise process as single", 2, 2, null, 1, 1).addFootnote(footId, footText1));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('c', "only-count", "Do not output reads"));
        optSet.addOpt(new Opt('t', "splitter-threads", "Number of splitter threads. No point setting too high, "
            + "i/o is the likely bottleneck and a writing thread will be spawned per each sample", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('o', "out-dir", "Output directory", 1).setDefaultValue("out_split"));
        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
        optSet.addOpt(new Opt(null, "append", "If output file(s) exist(s) for a given sample, append"));
        optSet.addOpt(new Opt(null, "force", "If output file(s) exist(s) for a given sample, force overwrite"));
        optSet.addOpt(new Opt('M', "[TODO] matchless-output", "Output reads with unmatched barcodes to R1/R2/SE file(s) prefixed with <arg>. If not set, these reads will be discarded", 1));
        footId++;
        String footText2 = "Consider increasing to sacrifice memory for speed. Decrease if encountering 'out of memory' errors.";
        optSet.addOpt(new Opt('u', "out-buffer-size", "Number of FASTQ records (reads or pairs) "
            + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText2));
        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
            2, 1, 32).addFootnote(footId, footText2));

        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void splitFiles(ArrayList<String> inputFilenamesList, OptSet optSet) {
        
        String restrictionRemnantI = (String) optSet.getOpt("restriction-remnant-I").getValueOrDefault();
        KeyMap keyMap = new KeyMap(KEY_FILE_NAME, TOOL_NAME, BLANK_SAMPLE_NAME, OUT_Q_CAPACITY, restrictionRemnantI);

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
        int k = -1;
        //READ INPUT AND POPULATE PairMers MAP
//            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);

        BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
//            boolean stranded = false;
        int ioThreads = keyMap.getSamplesTotal() + 1;
        ArrayList<Future<?>> ioFutures = new ArrayList<>(ioThreads);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(ioThreads, ioThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
//        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFilenamesList, k, TOOL_NAME, "records", IN_BUFFER_SIZE);
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFilenamesList, TOOL_NAME, "records", IN_BUFFER_SIZE);
        ioFutures.add(ioExecutorService.submit(inputReaderProducer));

        long timeStart = System.currentTimeMillis();
        int count = 0;
        while (inputReaderProducer.getGuessedInputFormat() == null) {
            try {
                //IF nothing happens after 5 seconds
                if (System.currentTimeMillis() - timeStart > 2500 && (count++ % 25 == 0)) {
                    Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", TOOL_NAME);
                }
                Thread.sleep(100); //wait for 1/10 of a second
            } catch (InterruptedException ex) {
            }
        }
//        //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
//        if (!inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.InFormat.KMERS)) {
//            Reporter.report("[FATAL]", "Only k-mer sets accepted as input. Guessed format: " + inputReaderProducer.getGuessedInputFormat(), getClass().getSimpleName());
//        } else {
        //Start KmergerConsumerProducer and OutputWriterConsumer threads
        final ExecutorService splitterExecutorService = new ThreadPoolExecutor(SPLITTER_THREADS, SPLITTER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> splittersFutures = new ArrayList<>(SPLITTER_THREADS);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(SPLITTER_THREADS * 5);
        for (int i = 0; i < SPLITTER_THREADS; i++) {
            splittersFutures.add(splitterExecutorService.submit(new SplitterConsumerProducer(inputQueue, keyMap, TOOL_NAME, OUT_BUFFER_SIZE,
                optSet, finalMessages)));
        }

        //WRITER THREADS not set if just counting
        if (!optSet.getOpt("only-count").getOptFlag()) {
            for (Map.Entry<String, BlockingQueue<PerSampleBuffer>> entrySet : keyMap.getSampleToQueueMap().entrySet()) {
                String sample = entrySet.getKey();
//                System.err.println("Starting writer thread for sample "+sample);
                BlockingQueue<PerSampleBuffer> queue = entrySet.getValue();
//            ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(OUT_DIR + "/" + sample, queue, TOOL_NAME, SPLITTER_THREADS)));
                ioFutures.add(ioExecutorService.submit(new FileWriterConsumer(queue, TOOL_NAME, OUT_DIR, sample, SPLITTER_THREADS, R1_SUFFIX, R2_SUFFIX, SE_SUFFIX, OVERWRITE, APPEND)));
            }
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
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
        for (Message fm : finalMessages) {
            Reporter.report(fm.getLevel().toString(), fm.getBody(), fm.getCaller());
        }
        if (optSet.getOpt("only-count").getOptFlag()) {
            ConcurrentHashMap<String, Long> sampleToCountMap = keyMap.getSampleToCountMap();
            Set<String> keySet = sampleToCountMap.keySet();
            for (String key : keySet) {
                Long c = sampleToCountMap.get(key);
                Reporter.report("[INFO]", NumberFormat.getInstance().format(c) + " records@sample " + key, TOOL_NAME);
            }
        }

    }

}
