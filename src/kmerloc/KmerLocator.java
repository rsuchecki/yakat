/*
 * Copyright 2020 rad.
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
package kmerloc;

import allorfs.AllOrfs;
import allorfs.OrfPredictorConsumer;
import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.io.PrintStream;
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
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import kmermatch.Kmer;
import kmermatch.KmerSetPopulatorConsumer;
import shared.FastaIndexed;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;
import shared.Sequence;
import shared.SequenceOps;
import shared.StdRedirect;

/**
 *
 * @author rad
 */
public class KmerLocator {

    //Read and populate k-mer set //later TODO: multiple k-mer sets - can process range of k values in single pass of the genome
    //Split target FASTA into separate sequences //Later possibly split sequences further
    //  
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 180;
    private final int IN_Q_CAPACITY;

    public KmerLocator(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        new StdRedirect(optSet, TOOL_NAME, StdRedirect.RedirectType.REDIRECT_OUT);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();

        locateKmers(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Identify and extract all (longest) ORFs from a genome");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('d', "domtbl", "The domtbl output generated by hmmer search ", 1).setRequired(true));
        optSet.addOpt(new Opt('f', "fasta", "Target FASTA file (indexed)", 1).setRequired(true));
        optSet.addOpt(new Opt('b', "first-bases", "Only consider first <arg> bases of a FASTA", 1).setMinValue(6));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
                2, 1, 256));
        int footId = 1;
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
//        optSet.addOpt(new Opt('v', "invert-matching", "Output unmatched reads"));
        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()));
        //.addFootnote(footId, "Caution! Large sequences times high number of threads => huge memory consumption"));

        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
//        optSet.addOpt(new Opt('s', "strand", "Output ORFs identified on strand ", 1).setDefaultValue("both"));
////        footId++;
        optSet.addOpt(new Opt('o', "out", "Print output to <arg> file", 1).setDefaultValue("/dev/stdout")); //on Windows use "CON" as default file name
        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_KMERS_FILENAMEs", "names of input files (each containing k-mers, optional \tcounts)", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void locateKmers(ArrayList<String> inputFilenamesList, OptSet optSet) {
        if (inputFilenamesList == null || inputFilenamesList.isEmpty()) { //READSTDIN
            inputFilenamesList = new ArrayList<>();
            inputFilenamesList.add("-");
        }
        //PARSE FASTA FAI INDEX
        String fasta = (String) optSet.getOpt("f").getValueOrDefault();
        FastaIndexed fastaIndexed = new FastaIndexed(TOOL_NAME, fasta, fasta + ".fai");

        int MAX_THREADS = (int) optSet.getOpt("t").getValueOrDefault();
        PrintStream bufferedOut = new PrintStream(new java.io.BufferedOutputStream(System.out, 65535));

        Integer firsBasesNum = null;
        if (optSet.getOpt("b").isUsed()) {
            firsBasesNum = (int) optSet.getOpt("b").getValueOrDefault();
        }

        BlockingQueue inputKmersQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        //SPAWN KMERS INPUT READING THREAD        
        ArrayList<Future<?>> inputKmersFutures = new ArrayList<>(1);
        final ExecutorService inputKmersExecutorService = new ThreadPoolExecutor(1, 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputKmersQueue, new ArrayList<Integer>(), inputFilenamesList, 8092, TOOL_NAME);

        inputKmersFutures.add(inputKmersExecutorService.submit(inputReaderProducer));

        Reporter.report("[INFO]", "Start populating k-mers set", TOOL_NAME);

        int THREADS = 2;
        Integer k = null;
        boolean storeASCII = true;

        //SPAWN MAP - POPULATOR THREADS
        ConcurrentSkipListSet<Kmer> kmers = new ConcurrentSkipListSet();
        final ExecutorService populatorExecutorService = new ThreadPoolExecutor(THREADS, THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> populatorFutures = new ArrayList<>(THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(THREADS * 5);
        for (int i = 0; i < THREADS; i++) {
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

        Pattern nonNuclPattern = Pattern.compile(".*[^acgtACGT]+.*");
        Pattern tab = Pattern.compile("\t");
        for (String id : fastaIndexed.getIds()) {
            Reporter.report("[INFO]", "Reading in " + id, TOOL_NAME);
            Sequence sequence = firsBasesNum == null ? new Sequence(id, fastaIndexed.getSequence(id)) : new Sequence(id, fastaIndexed.getSequence(id, 1L, firsBasesNum.longValue()));

//            System.out.println(sequence.getLength());
            //KMERIZE
            
            int KEY = inputReaderProducer.getKmerLengths().get(0);
            
            int maxKmer = sequence.getLength() - KEY + 1;
            int count = 0;
            for (int i = 0; i < maxKmer; i++) {
                String canonical = SequenceOps.getCanonical(sequence.getSequenceString().subSequence(i, i + KEY ).toString());
                if (kmers.contains(new Kmer(canonical, storeASCII))) {
                    count++;
                    bufferedOut.append(id+":"+(i+1)+"-"+(i+KEY)+"\t"+canonical).append(System.lineSeparator());
                }
            }
            Reporter.report("[INFO]", "Matched "+count+" kmers in " + id, TOOL_NAME);
        }

//                
//        final ExecutorService splitterExecutorService = new ThreadPoolExecutor(MAX_THREADS, MAX_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
//        ArrayList<Future<?>> splittersFutures = new ArrayList<>(MAX_THREADS);
////        AtomicInteger splitterThreads = new AtomicInteger(MAX_THREADS);
//        ArrayList<Message> finalMessages = new ArrayList<>(MAX_THREADS * 5);
//        BlockingQueue<ArrayList<Sequence>> inputQueue = new ArrayBlockingQueue<>(2);
//        for (int i = 0; i < MAX_THREADS; i++) {
//            splittersFutures.add(splitterExecutorService.submit(new OrfPredictorConsumer(inputQueue, TOOL_NAME, minLength, bufferedOut, 
//                    translate, requireStop, strand)));
//        }
////        try (PrintStream bufferedOut = new PrintStream(new java.io.BufferedOutputStream(System.out, 65535))) {
//        try {
//            for (String id : fastaIndexed.getIds()) {
//                Reporter.report("[INFO]", "Reading in " + id, TOOL_NAME);
////                Sequence sequence = new Sequence(id, fastaIndexed.getSequence(id, 1L, 50000000L));
//                Sequence sequence = firsBasesNum == null ? new Sequence(id, fastaIndexed.getSequence(id)) : new Sequence(id, fastaIndexed.getSequence(id, 1L, firsBasesNum.longValue()));
//                ArrayList<Sequence> seqList = new ArrayList<>();
//                seqList.add(sequence);
//                inputQueue.put(seqList);
//                //                Reporter.report("[INFO]", "Identifying ORFs in " + id, TOOL_NAME);
//                //                ArrayList<Orf> orfs = sequence.getOrfs(minLength);
//                //                Reporter.report("[INFO]",  NumberFormat.getInstance().format(orfs.size())+" ORFs found ", TOOL_NAME);
//                //                for (Orf orf : orfs) {
//                ////                System.out.printf("%8s%12d%12d%3d%12d\n",orf.getParenId(),orf.getFrom(),orf.getTo(),orf.getFrame(),orf.getLength());
//                //                    if (orf.hasStopCodon()) {
//                //                        CharSequence seq = sequence.getSequenceString().subSequence(orf.getFrom() - 1, orf.getTo());
//                //                        bufferedOut.append(orf.getFastaHeader());
//                //                        bufferedOut.append(System.lineSeparator());
//                //                        if (orf.getFrame() < 0) {
//                //                            seq = SequenceOps.getReverseComplement(seq);
//                //                        }
//                //                        bufferedOut.append(SequenceOps.translate(seq));
//                //                        bufferedOut.append(System.lineSeparator());
//                //                    }
//                //                }
//            }
//            inputQueue.put(new ArrayList<>());
//        } catch (InterruptedException ex) {
//            Logger.getLogger(AllOrfs.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//        splitterExecutorService.shutdown();
//        try {
//            for (Future<?> f : splittersFutures) {
//                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            }
//            splitterExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
//        } catch (ExecutionException ex) {
//            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
//            ex.printStackTrace();
//        } catch (TimeoutException ex) {
//            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
//        }
        bufferedOut.flush();

        Reporter.report("[INFO]", "Finished extracting ORFs", TOOL_NAME);
    }
}
