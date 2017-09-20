/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kmermatch;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerMatch {

    private final String TOOL_NAME;
    private int MATCHER_THREADS = 1;

    //IN
    private final int IN_BUFFER_SIZE;
    private final int IN_Q_CAPACITY;

    //OUT
    private final int OUT_BUFFER_SIZE;
    private final int OUT_Q_CAPACITY;
    private final int HELP_WIDTH = 180;

    public KmerMatch(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();

        MATCHER_THREADS = (int) optSet.getOpt("t").getValueOrDefault();

        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
        OUT_Q_CAPACITY = (int) optSet.getOpt("q").getValueOrDefault();

        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        match(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Given a set of k-mers and (by default one-per-line) FASTQ input, output FASTQ records if sufficient matching k-mers are present ");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('K', "k-mers", "Set of k-mers to be used as reference for matching the reads", 1).setRequired(true));
//        optSet.addOpt(new Opt('A', "expected-adapter", "Expected adapter sequence (a fragment will do)", 1).setDefaultValue("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"));
//        optSet.addOpt(new Opt(null, "adapter-prefix-length", "Length of the adapter prefix used to identify 3' read-through", 1).setDefaultValue(9));
//        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to in-queue", 1024, 1, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
                2, 1, 256));
//        //TRIMMING AND LENGTH
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Trimming and length settings]");
//        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
//        optSet.addOpt(new Opt('a', "keep-adapters", "Do not trim adapters found next to PstI and MspI sites "));
//        optSet.addOpt(new Opt('p', "keep-non-PstI-starting", "Keep reads (or pairs) which do not start with 'barcodeTGCAG'"));
//
        int footId = 1;
//        String footText1 = "Note that certain combinations of min-length-* settings can lead to both mates of a pair ending up in SE/orphans output file.";
//        optSet.addOpt(new Opt('r', "min-length-single-read", "Only output a read if length is no less than <arg> bp", 1, 1, null, 1, 1).addFootnote(footId, footText1));
//        optSet.addOpt(new Opt('e', "min-length-pair-each", "Only output a read pair if length of each is no less than <arg> bp, otherwise process as single", 1, 1, null, 1, 1).addFootnote(footId, footText1));
//        optSet.addOpt(new Opt('s', "min-length-pair-sum", "Only output a read pair if combined length is no less than <arg> bp, otherwise process as single", 2, 2, null, 1, 1).addFootnote(footId, footText1));
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('v', "invert-matching", "Output unmatched reads"));
        optSet.addOpt(new Opt('t', "threads", "Number of threads. ", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
//        optSet.addOpt(new Opt('o', "out-dir", "Output directory", 1).setDefaultValue("out_split"));
//        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
//        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
//        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
//        optSet.addOpt(new Opt(null, "append", "If output file(s) exist(s) for a given sample, append"));
//        optSet.addOpt(new Opt(null, "force", "If output file(s) exist(s) for a given sample, force overwrite"));
//        optSet.addOpt(new Opt('M', "[TODO] matchless-output", "Output reads with unmatched barcodes to R1/R2/SE file(s) prefixed with <arg>. If not set, these reads will be discarded", 1));
//        footId++;
        optSet.addOpt(new Opt('o', "out-file", "Send output to <arg> file", 1).setDefaultValue("/dev/stdout"));

        String footText2 = "Consider increasing to sacrifice memory for speed. Decrease if encountering 'out of memory' errors.";
        optSet.addOpt(new Opt('u', "out-buffer-size", "Number of FASTQ records (reads or pairs) "
                + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText2));
        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
                2, 1, 256).addFootnote(footId, footText2));

        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void match(ArrayList<String> inputFilenamesList, OptSet optSet) {
        ArrayList<String> inputFilesToMatch = new ArrayList<>();
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilesToMatch.addAll(po.getValues());
            }
        }
        //READ INPUT AND POPULATE IDs MAP
        
        
        
        BlockingQueue inputKmersQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);

        //SPAWN IDS INPUT READING THREAD        
        ArrayList<Future<?>> inputKmersFutures = new ArrayList<>(1);
        final ExecutorService inputKmersExecutorService = new ThreadPoolExecutor(1, 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        InputReaderProducer inputReaderProducer  = new InputReaderProducer(inputKmersQueue, new ArrayList<Integer>(), optSet.getOpt("-K").getValues() , IN_BUFFER_SIZE, TOOL_NAME); 
        
        inputKmersFutures.add(inputKmersExecutorService.submit(inputReaderProducer));

        

        //SPAWN MAP - POPULATOR THREADS
        ConcurrentSkipListSet<Kmer> kmers = new ConcurrentSkipListSet();
        final ExecutorService populatorExecutorService = new ThreadPoolExecutor(MATCHER_THREADS, MATCHER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> populatorFutures = new ArrayList<>(MATCHER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(MATCHER_THREADS * 5);
        for (int i = 0; i < MATCHER_THREADS; i++) {
            populatorFutures.add(populatorExecutorService.submit(new KmerSetPopulatorConsumer(kmers, inputKmersQueue)));
        }
        populatorExecutorService.shutdown();
        inputKmersExecutorService.shutdown();
        try {

            for (Future<?> f : populatorFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            populatorExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            for (Future<?> f : inputKmersFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            inputKmersExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }

        Reporter.report("[INFO]", "Finished populating kmers set, n=" + NumberFormat.getNumberInstance().format(kmers.size()), TOOL_NAME);

        
             
        
        
        //NOW PROCESS INPUT READS
        String outFile = (String) optSet.getOpt("out-file").getValueOrDefault();
        BlockingQueue outputQueue = new ArrayBlockingQueue(OUT_Q_CAPACITY);

        BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        ArrayList<Future<?>> ioFutures = new ArrayList<>(2);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(2, 2, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        //READER THREAD
        InputReaderProducer inputReaderProducer2 = new InputReaderProducer(inputQueue, inputFilenamesList, TOOL_NAME, "records", IN_BUFFER_SIZE);
        ioFutures.add(ioExecutorService.submit(inputReaderProducer2));

        
        long timeStart = System.currentTimeMillis();
        int count = 0;
        while (inputReaderProducer2.getGuessedInputFormat() == null) {
            try {
                //IF nothing happens after 5 seconds
                if (System.currentTimeMillis() - timeStart > 2500 && (count++ % 25 == 0)) {
                    Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", getClass().getSimpleName());
                }
                Thread.sleep(100); //wait for 1/10 of a second
            } catch (InterruptedException ex) {
            }
        }
        
        
        
        //WRITER THREAD
        ioFutures.add(ioExecutorService.submit(new WriterConsumer(outputQueue, outFile, MATCHER_THREADS, TOOL_NAME)));
//        ioFutures.add(ioExecutorService.submit(new shared.WriterConsumer(outputQueue, "MATCHED_FASTQ", TOOL_NAME)));

        //SPAWN MATCHER-THREADS
        final ExecutorService matcherExecutorService = new ThreadPoolExecutor(MATCHER_THREADS, MATCHER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> matcherFutures = new ArrayList<>(MATCHER_THREADS);
        for (int i = 0; i < MATCHER_THREADS; i++) {            
            matcherFutures.add(matcherExecutorService.submit(new KmerMatcherConsumerProducer(inputQueue, outputQueue, kmers, IN_BUFFER_SIZE, TOOL_NAME, 
                    finalMessages, optSet.getOpt("v").isUsed(), inputReaderProducer.getGuessedInputFormat(), inputReaderProducer.getKmerLengths().get(0))));
        }

        matcherExecutorService.shutdown();
        ioExecutorService.shutdown();

        try {

            for (Future<?> f : matcherFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            matcherExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
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

    }

}
