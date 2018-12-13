/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fastqmatchid;

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
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class FastqMatchId {

    private final String TOOL_NAME;
    private int MATCHER_THREADS = 1;

    //IN
    private final int IN_BUFFER_SIZE;
    private final int IN_Q_CAPACITY;

    //OUT
    private final int OUT_BUFFER_SIZE;
    private final int OUT_Q_CAPACITY;
    private final int HELP_WIDTH = 180;

    public FastqMatchId(String[] args, String callerName, String toolName) {
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
        matchIds(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Given a set of read identifiers and (by default one-per-line) FASTQ input, output FASTQ records for the matching ids");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('I', "identifiers", "Read identifiers to be matched ", 1).setRequired(true));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to in-queue", 1024, 1, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
                2, 1, 256));
        int footId = 1;
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('v', "invert-matching", "Output unmatched reads"));
        optSet.addOpt(new Opt('t', "threads", "Number of threads. ", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");

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

    private void matchIds(ArrayList<String> inputFilenamesList, OptSet optSet) {
        //READ INPUT AND POPULATE IDs MAP

        BlockingQueue inputIdsQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);

        //SPAWN IDS INPUT READING THREAD        
        ArrayList<Future<?>> inputIdsFutures = new ArrayList<>(1);
        final ExecutorService inputIdsExecutorService = new ThreadPoolExecutor(1, 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputIdsQueue, (String) optSet.getOpt("-I").getValueIfSingle(), TOOL_NAME, "records", IN_BUFFER_SIZE, true);
        inputIdsFutures.add(inputIdsExecutorService.submit(inputReaderProducer));

        

        //SPAWN MAP - POPULATOR THREADS
        ConcurrentSkipListSet<Identifier> ids = new ConcurrentSkipListSet();
        final ExecutorService populatorExecutorService = new ThreadPoolExecutor(MATCHER_THREADS, MATCHER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> populatorFutures = new ArrayList<>(MATCHER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(MATCHER_THREADS * 5);
        for (int i = 0; i < MATCHER_THREADS; i++) {
            populatorFutures.add(populatorExecutorService.submit(new IdSetPopulatorConsumer(ids, inputIdsQueue)));
        }
        populatorExecutorService.shutdown();
        inputIdsExecutorService.shutdown();
        try {

            for (Future<?> f : populatorFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            populatorExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            for (Future<?> f : inputIdsFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            inputIdsExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }

        Reporter.report("[INFO]", "Finished populating identifiers set, n=" + NumberFormat.getNumberInstance().format(ids.size()), TOOL_NAME);

        
             
        
        
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
                    Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", TOOL_NAME);
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
            matcherFutures.add(matcherExecutorService.submit(new MatcherConsumerProducer(inputQueue, outputQueue, ids, IN_BUFFER_SIZE, TOOL_NAME, finalMessages, optSet.getOpt("v").isUsed())));
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
