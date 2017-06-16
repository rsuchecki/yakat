/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hmmerdoms;

import fastqmatchid.*;
import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import gbssplit.Sample;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import shared.FastaIndexed;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class HmmerDomains {

    private final String TOOL_NAME;
//    private int MATCHER_THREADS = 1;

    //IN
    private final int IN_BUFFER_SIZE;
    private final int IN_Q_CAPACITY;

    
    //OUT
//    private final int OUT_BUFFER_SIZE;
//    private final int OUT_Q_CAPACITY;
    private final int HELP_WIDTH = 180;

    public HmmerDomains(String[] args, String callerName, String toolName) {
        ArrayList<String> inputFilenamesList = new ArrayList<>();
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        TOOL_NAME = callerName + " " + toolName;
        //PARSE OPTS
        IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();
        
        
//        MATCHER_THREADS = (int) optSet.getOpt("t").getValueOrDefault();
//        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
//        OUT_Q_CAPACITY = (int) optSet.getOpt("q").getValueOrDefault();
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFilenamesList.addAll(po.getValues());
            }
        }
        processHmmerDomains(inputFilenamesList, optSet);

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Given tabular .domtbl output from hmmmer search against six-frame translation of a genome, extract putative genomic regions coordinates."
                + "\n(1) Second and third \"_\" delimited field in identifier reflect the offset introduced when conigerizing a chromosome and the number of removed Ns respectively"
                + "\n(2) Expecting translation frame in FASTA description \\in {fr1,fr2,...,fr6} (following FASTA id and whitespace)");
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('d', "domtbl", "The domtbl output generated by hmmer search ", 1).setRequired(true));
        optSet.addOpt(new Opt('f', "fasta", "Reference FASTA file ", 1).setRequired(true));
//        optSet.addOpt(new Opt('A', "expected-adapter", "Expected adapter sequence (a fragment will do)", 1).setDefaultValue("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"));
//        optSet.addOpt(new Opt(null, "adapter-prefix-length", "Length of the adapter prefix used to identify 3' read-through", 1).setDefaultValue(9));
//        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
                + "passed to in-queue", 1024, 128, 8092));
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
//        optSet.addOpt(new Opt('v', "invert-matching", "Output unmatched reads"));
//        optSet.addOpt(new Opt('t', "threads", "Number of splitter threads. No point setting too high, "
//                + "i/o is the likely bottleneck and a writing thread will be spawned per each sample", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));

        //OUTPUT
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
        optSet.addOpt(new Opt('g', "max-gap", "Maximum bp gap allowed between grouped domain-hits", 1).setDefaultValue(200));
//        optSet.addOpt(new Opt('e', "min-elems", "Minimum elements for a group of domain-hits to be reported", 1).setDefaultValue(1).setMinValue(1));
        
//        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
//        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
//        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
//        optSet.addOpt(new Opt(null, "append", "If output file(s) exist(s) for a given sample, append"));
//        optSet.addOpt(new Opt(null, "force", "If output file(s) exist(s) for a given sample, force overwrite"));
//        optSet.addOpt(new Opt('M', "[TODO] matchless-output", "Output reads with unmatched barcodes to R1/R2/SE file(s) prefixed with <arg>. If not set, these reads will be discarded", 1));
//        footId++;
        optSet.addOpt(new Opt('o', "out-file", "Print output to <arg> file", 1).setDefaultValue("/dev/stdout")); //on Windows use "CON" as default file name

//        String footText2 = "Consider increasing to sacrifice memory for speed. Decrease if encountering 'out of memory' errors.";
//        optSet.addOpt(new Opt('u', "out-buffer-size", "Number of records  "
//                + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText2));
//        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
//                2, 1, 256).addFootnote(footId, footText2));
        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files (hmmer search domtbl)", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

    
    private void processHmmerDomains(ArrayList<String> inputFilenamesList, OptSet optSet) {

        
        DomainHitsPerTarget domainHitsPerTarget = new DomainHitsPerTarget();
        
        
        String inputFile = inputFilenamesList.get(0);

        Pattern spliPattern = Pattern.compile("\t| +");
        BufferedReader content = null;
        try {
            if (inputFile.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile));
                content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"));
            } else {
                content = new BufferedReader(new FileReader(new File(inputFile)));
            }
            String line;
            while ((line = content.readLine()) != null && !line.isEmpty()) {
//                String[] toks = line.split("\t");
                if (!line.startsWith("#")) {
                    String[] toks = spliPattern.split(line);
                    DomainHit domainHit = new DomainHit(toks);
                    domainHitsPerTarget.addHit(domainHit);
//                    System.out.println(line);
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found exception: " + ex.getMessage(), TOOL_NAME);
            System.exit(1);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        } finally {
            try {
                if (content != null) {
                    content.close();
                }
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
        FastaIndexed fastaIndexed = null;
        if(optSet.getOpt("f").isUsed()) {
            String fasta = (String) optSet.getOpt("f").getValueOrDefault();
            fastaIndexed = new FastaIndexed(TOOL_NAME, fasta, fasta+".fai");
        }
//        System.err.println(fastaIndexed.getSequence("LOC_Os10g35790.1", null, null));
        
        domainHitsPerTarget.processDomainsPerStrand(true,(int)optSet.getOpt("g").getValueOrDefault(), fastaIndexed);
        domainHitsPerTarget.processDomainsPerStrand(false,(int)optSet.getOpt("g").getValueOrDefault(), fastaIndexed);

    }

    
 }
