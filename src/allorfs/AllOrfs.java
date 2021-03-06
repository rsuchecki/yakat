/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package allorfs;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.io.PrintStream;
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
import shared.FastaIndexed;
import shared.Message;
import shared.Reporter;
import shared.Sequence;
import shared.StdRedirect;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class AllOrfs {

    private final String TOOL_NAME;
//    private int MATCHER_THREADS = 1;

    //IN
//    private final int IN_BUFFER_SIZE;
//    private final int IN_Q_CAPACITY;

    //OUT
    private final int WRITER_BUFFER_SIZE = 1024;
//    private final int OUT_BUFFER_SIZE;
//    private final int OUT_Q_CAPACITY;
    private final int HELP_WIDTH = 180;

    public AllOrfs(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
//        IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
//        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();

//        MATCHER_THREADS = (int) optSet.getOpt("t").getValueOrDefault();
//        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
//        OUT_Q_CAPACITY = (int) optSet.getOpt("q").getValueOrDefault();
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


        allOrfs(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Identify and extract all (longest) ORFs from a genome");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('d', "domtbl", "The domtbl output generated by hmmer search ", 1).setRequired(true));
        optSet.addOpt(new Opt('f', "fasta", "Reference FASTA file ", 1).setRequired(true));
        optSet.addOpt(new Opt('n', "no-strict-stop-codon", "Output putative ORFs even when missing a stop codon"));
        optSet.addOpt(new Opt('d', "do-not-translate", "Output putative ORFs as  nucleotide sequences"));

//        //TRIMMING AND LENGTH
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Trimming and length settings]");
        optSet.addOpt(new Opt('b', "first-bases", "Only consider first <arg> bases of a FASTA entry for ORF prediction", 1).setMinValue(6));
        int footId = 1;
        //RUNTIME
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
//        optSet.addOpt(new Opt('v', "invert-matching", "Output unmatched reads"));
        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(footId, "Caution! Large sequences times high number of threads => huge memory consumption"));

        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('m', "min-orf-length", "Minimum output ORF length", 1).setDefaultValue(100).setMinValue(6));
        optSet.addOpt(new Opt('s', "strand", "Output ORFs identified on strand ", 1).setDefaultValue("both"));
////        footId++;
        optSet.addOpt(new Opt('o', "out", "Print output to <arg> file", 1).setDefaultValue("/dev/stdout")); //on Windows use "CON" as default file name
        //POSITIONAL
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files (hmmer search domtbl)", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    private void allOrfs(ArrayList<String> inputFilenamesList, OptSet optSet) {
        //PARSE INPUT FROm HMMER
//        DomainHitsPerTarget domainHitsPerTarget = new DomainHitsPerTarget();
//        Pattern spliPattern = Pattern.compile("\t| +");
//        BufferedReader content = null;

        if (inputFilenamesList == null || inputFilenamesList.isEmpty()) { //READSTDIN
            inputFilenamesList = new ArrayList<>();
            inputFilenamesList.add("-");
        }
        //PARSE FASTA FAI INDEX
        String fasta = (String) optSet.getOpt("f").getValueOrDefault();
        FastaIndexed fastaIndexed = new FastaIndexed(TOOL_NAME, fasta, fasta + ".fai");

//        boolean greedyExtendIfIUnterminatedORF = true;
//        if (optSet.getOpt("n").isUsed()) {
//            greedyExtendIfIUnterminatedORF = false;
//        }
        Integer firsBasesNum = null;
        if (optSet.getOpt("b").isUsed()) {
            firsBasesNum = (int) optSet.getOpt("b").getValueOrDefault();
        }
        int minLength = (int) optSet.getOpt("m").getValueOrDefault();
        boolean requireStop = !optSet.getOpt("n").isUsed();
        boolean translate = !optSet.getOpt("d").isUsed();
        String strand =  (String) optSet.getOpt("s").getValueOrDefault();
        int MAX_THREADS = (int) optSet.getOpt("t").getValueOrDefault();
        PrintStream bufferedOut = new PrintStream(new java.io.BufferedOutputStream(System.out, 65535));

        final ExecutorService splitterExecutorService = new ThreadPoolExecutor(MAX_THREADS, MAX_THREADS, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> splittersFutures = new ArrayList<>(MAX_THREADS);
//        AtomicInteger splitterThreads = new AtomicInteger(MAX_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(MAX_THREADS * 5);
        BlockingQueue<ArrayList<Sequence>> inputQueue = new ArrayBlockingQueue<>(2);
        for (int i = 0; i < MAX_THREADS; i++) {
            splittersFutures.add(splitterExecutorService.submit(new OrfPredictorConsumer(inputQueue, TOOL_NAME, minLength, bufferedOut, 
                    translate, requireStop, strand)));
        }
//        try (PrintStream bufferedOut = new PrintStream(new java.io.BufferedOutputStream(System.out, 65535))) {
        try {
            for (String id : fastaIndexed.getIds()) {
                Reporter.report("[INFO]", "Reading in " + id, TOOL_NAME);
//                Sequence sequence = new Sequence(id, fastaIndexed.getSequence(id, 1L, 50000000L));
                Sequence sequence = firsBasesNum == null ? new Sequence(id, fastaIndexed.getSequence(id)) : new Sequence(id, fastaIndexed.getSequence(id, 1, firsBasesNum));
                ArrayList<Sequence> seqList = new ArrayList<>();
                seqList.add(sequence);
                inputQueue.put(seqList);
                //                Reporter.report("[INFO]", "Identifying ORFs in " + id, TOOL_NAME);
                //                ArrayList<Orf> orfs = sequence.getOrfs(minLength);
                //                Reporter.report("[INFO]",  NumberFormat.getInstance().format(orfs.size())+" ORFs found ", TOOL_NAME);
                //                for (Orf orf : orfs) {
                ////                System.out.printf("%8s%12d%12d%3d%12d\n",orf.getParenId(),orf.getFrom(),orf.getTo(),orf.getFrame(),orf.getLength());
                //                    if (orf.hasStopCodon()) {
                //                        CharSequence seq = sequence.getSequenceString().subSequence(orf.getFrom() - 1, orf.getTo());
                //                        bufferedOut.append(orf.getFastaHeader());
                //                        bufferedOut.append(System.lineSeparator());
                //                        if (orf.getFrame() < 0) {
                //                            seq = SequenceOps.getReverseComplement(seq);
                //                        }
                //                        bufferedOut.append(SequenceOps.translate(seq));
                //                        bufferedOut.append(System.lineSeparator());
                //                    }
                //                }
            }
            inputQueue.put(new ArrayList<>());
        } catch (InterruptedException ex) {
            Logger.getLogger(AllOrfs.class.getName()).log(Level.SEVERE, null, ex);
        }

        splitterExecutorService.shutdown();
        try {
            for (Future<?> f : splittersFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            splitterExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
        bufferedOut.close();

        Reporter.report("[INFO]", "Finished extracting ORFs", TOOL_NAME);
    }
//    }

}
