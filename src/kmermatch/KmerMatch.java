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
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
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
import kmerextender.CoreCoder;
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
        OptSet optSet = new OptSet("Given a set of k-mers (or FASTA/FASTQ to be k-merized) and (by default one-per-line) FAST[A|Q] input, output FAST[A|Q] records if sufficient matching k-mers are present. Do not mix different input types in one run (FASTA/FASTQ/PE/SE) as input format is recognized early on and subsequent records are assumed to be in the same format.");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('K', "k-mers", "Sequence(s) to be used as reference for matching the reads, by default a set of k-mers is expected. If <arg> needs to be k-merized then -k / --k-mer-length is also required.", 1).setRequired(true));
        optSet.addOpt(new Opt('k', "k-mer-length", "Specify required k-mer size if reference is to be k-merized", 1).setMinValue(2).setMaxValue(1024));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to in-queue", 2048, 1, 8092));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
                2, 1, 256));
        int footId = 1;
//        String footText1 = "Note that certain combinations of min-length-* settings can lead to both mates of a pair ending up in SE/orphans output file.";
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
        optSet.addOpt(new Opt('m', "min-matches", "Minimum number of k-mers from the targeted sequence (pair) present in the reference",1).setMinValue(0).setDefaultValue(1));
        optSet.addOpt(new Opt('M', "min-matches-fraction", "Minimum fraction of k-mers from the targeted sequence (pair) present in the reference", 0.0, 0.0, 1.0));
        optSet.addOpt(new Opt('v', "invert-matching", "Invert matching: output targets that contain fewer than -m <arg> k-mers matching the reference set"));
        optSet.addOpt(new Opt('a', "ascii-encoding", "Store k-mers in plain ascii - possibly faster than the default encoding, but will consume more memory."));
        
//                optSet.addOpt(new Opt('S', "stranded-matching", "Do not reverse-complement k-mers (by default canonical representation of a k-mer is stored and matched)"));
//        optSet.addOpt(new Opt('B', "match-both-mates", "Relevant for PE input only. Demand each mate to have -m <arg> k-mers matching the reference set, by default, both mates are caught if at least one has [-M] kmer(s) matching the reference"));
//        optSet.addOpt(new Opt('A', "report-all", "Report all input sequences, with the number of matching k-mers appended to the identifier, other matching settings will be ignored"));

        
        
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
        optSet.addOpt(new Opt('d', "dump-kmers", "Dump reference k-mers to stdout - implemented for debugging purposes only so not optimized"));

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
        //READ INPUT AND POPULATE IDs MAP
        boolean storeASCII = optSet.getOpt("-a").isUsed();        
        Integer k = (Integer) optSet.getOpt("-k").getValueIfSingle();
        ArrayList<Integer> kValues = new ArrayList<>();
        kValues.add(k != null ? k : 0);
        
        BlockingQueue inputKmersQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        //SPAWN IDS INPUT READING THREAD        
        ArrayList<Future<?>> inputKmersFutures = new ArrayList<>(1);
        final ExecutorService inputKmersExecutorService = new ThreadPoolExecutor(1, 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        InputReaderProducer inputReaderProducer = new InputReaderProducer( inputKmersQueue, kValues, optSet.getOpt("-K").getValues(), IN_BUFFER_SIZE, TOOL_NAME);

        inputKmersFutures.add(inputKmersExecutorService.submit(inputReaderProducer));
        
        Reporter.report("[INFO]", "Start populating k-mers set", TOOL_NAME);

        //SPAWN MAP - POPULATOR THREADS
        ConcurrentSkipListSet<Kmer> kmers = new ConcurrentSkipListSet();
        final ExecutorService populatorExecutorService = new ThreadPoolExecutor(MATCHER_THREADS, MATCHER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> populatorFutures = new ArrayList<>(MATCHER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(MATCHER_THREADS * 5);
        for (int i = 0; i < MATCHER_THREADS; i++) {
            populatorFutures.add(populatorExecutorService.submit(new KmerSetPopulatorConsumer(kmers, inputKmersQueue, k, storeASCII)));
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

        Reporter.report("[INFO]", "Finished populating k-mers set, n=" + NumberFormat.getNumberInstance().format(kmers.size()), TOOL_NAME);

//        for(String inputFileName: inputFilenamesList) {
//            
//        }
        //NOW PROCESS INPUT READS, BUT FIRST SET-UP OUTPUT WRITING 
        String outFile = (String) optSet.getOpt("out-file").getValueOrDefault();
        BlockingQueue outputQueue = new ArrayBlockingQueue(OUT_Q_CAPACITY);
        BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        ArrayList<Future<?>> ioFutures = new ArrayList<>(2);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(2, 2, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        //READER THREAD
        InputReaderProducer inputReaderProducer2 = new InputReaderProducer(inputQueue, inputFilenamesList, TOOL_NAME, "records", IN_BUFFER_SIZE);
//        InputReaderProducer inputReaderProducer2 = new InputReaderProducer(inputQueue, kValues, inputFilenamesList, IN_BUFFER_SIZE, TOOL_NAME);
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

        k = k == null ? inputReaderProducer.getKmerLengths().get(0) : k;

        //WRITER THREAD
        ioFutures.add(ioExecutorService.submit(new WriterConsumer(outputQueue, outFile, MATCHER_THREADS, TOOL_NAME)));
//        ioFutures.add(ioExecutorService.submit(new shared.WriterConsumer(outputQueue, "MATCHED_FASTQ", TOOL_NAME)));

        //SPAWN MATCHER-THREADS
        final ExecutorService matcherExecutorService = new ThreadPoolExecutor(MATCHER_THREADS, MATCHER_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> matcherFutures = new ArrayList<>(MATCHER_THREADS);
        for (int i = 0; i < MATCHER_THREADS; i++) {
            matcherFutures.add(matcherExecutorService.submit(new KmerMatcherConsumerProducer(inputQueue, outputQueue, kmers, IN_BUFFER_SIZE, TOOL_NAME,
                    finalMessages, optSet.getOpt("v").isUsed(), (int) optSet.getOpt("m").getValueOrDefault(),
                    (double) optSet.getOpt("M").getValueOrDefault(),
                    inputReaderProducer.getGuessedInputFormat(), k, storeASCII)));
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
        
        
        if(optSet.getOpt("d").isUsed()) {
            Iterator<Kmer> iterator = kmers.iterator();
            while (iterator.hasNext()) {
                KmerBytes next = (KmerBytes) iterator.next();
                if(storeASCII) {
                    System.out.println(new String(next.getBytes(), StandardCharsets.US_ASCII));                                            
                } else {
                    System.out.println(CoreCoder.decodeCore(k, next.getBytes()));                                                                
                }
            }
        }
    }

}
