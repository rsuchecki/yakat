/*
 * Copyright 2016 rad.
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
package freqmers;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
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
public class FreqMers {

    private final static String DELIMITER = "\t";
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int REPORTING_SHIFT = 6; //input progress reporting every READER_BUFFER_SIZE*REPORTING_FACTOR records
    private final int READER_BUFFER_SIZE = 8192;
    private final int WRITER_BUFFER_SIZE = 8192;
    private HashMap<String, ArrayList<KmerLink>> map;
    private ArrayList<KmerFilter> kmerFilters;

    private boolean DEBUG = false;

//    private enum OutFmt {
//        IUPAC, AB, SLASH, FREQ, COV, EVIDENCE, CALL_WITH_EVIDENCE;
//    }

//    }
    public FreqMers(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        new StdRedirect(optSet, TOOL_NAME);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }
        if (optSet.getOpt("D").isUsed()) {
            DEBUG = true;
        }
        //READ IN PREVIOUSLY IDENTIFIED SET OF SNPs 
        //AND MAP SNP-OVERLAPPING k-mers TO SNP OBJECT
        buildFreqMerMap(optSet);
//        removeSnpMersPoorlyCoveredInParents(optSet); //DONE AFTER FALSE POSITIVE FILTERING

//        int minTotal = (int) optSet.getOpt("min-k-mer-frequency-sum").getValueOrDefault();
//        int minMinor = (int) optSet.getOpt("min-k-mer-frequency-minor").getValueOrDefault();
//        double minCoverage = (double) optSet.getOpt("min-snp-coverage").getValueOrDefault();
//        double maxError = (double) optSet.getOpt("max-coverage-error").getValueOrDefault();
        ArrayList<String> kmersFileNames = (ArrayList<String>) optSet.getOpt("K").getValues();
        outputFasta((String) optSet.getOpt("out-fasta").getValueOrDefault());
        ArrayList<String> sampleNames = threadKmersThroughMap(optSet, kmersFileNames);
        reportResults(sampleNames, (String) optSet.getOpt("out-stats").getValueOrDefault());
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("[A] Take FASTA sequences "
            + "[B] Take a set of k-mers per sample of interest. "
            + "[C] Record frequencies of k-mers overlapping the sequences "
            + "[D] Report k-mer coverage and frequencies along sequences");

        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('k', "k-mer-length", "It must match the size of the k-mers used to query the niks-snps", 1).setRequired(true).setMinValue(3).setMaxValue(255));
        optSet.addOpt(new Opt('f', null, "FASTA file", 1).setRequired(true));
//        optSet.addOpt(new Opt('f', "niks-fasta", "The (msa) FASTA file matching the SNP information", 1));
        optSet.addOpt(new Opt('K', "per-sample-k-mers", "A set (or sets) of k-mers to be threaded through the map of k-mer-links to SNPs").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 8192, 128, 65535));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up", 64, 1, 256));

        int footId = 1;
        String footText = "k-mer frequency corresponding to a SNP allele is obtained by taking a median of frequencies of all k-mers overlapping the underlying SNP";
        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Variant calling and reporting]");
//        optSet.addOpt(new Opt(null, "min-k-mer-frequency-sum", "Minimum frequency of k-mers which overlap "
//                + "with a SNP (minor+major allele)", 1).setMinValue(1).setDefaultValue(3).addFootnote(footId, footText));
//        optSet.addOpt(new Opt(null, "min-k-mer-frequency-minor", "Minimum frequency of k-mers "
//                + "which support the minor allele", 1).setMinValue(1).setDefaultValue(3).addFootnote(footId, footText));
////        optSet.addOpt(new Opt(null, "min-overlapping-k-mers", "At least <arg> k-mers must overlap a locus (for each allele)", 1).setMinValue(1).setDefaultValue(1));
//        optSet.addOpt(new Opt(null, "min-snp-coverage", "At least <arg> fraction of unique snp-covering k-mers must be present in a genotyped dataset (for each allele)", 1).setDefaultValue(0.9));
//        optSet.addOpt(new Opt(null, "max-coverage-error", "Maximum fraction of unique snp-covering k-mers that will be ignored as errorneous assignment and not considered for genotyping", 1).setMinValue(0.00).setDefaultValue(0.05));
//
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime and output settings]");
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
        optSet.addOpt(new Opt(null, "out-fasta", "Output annotated sequences to this file (in FASTA format)", 1));
        optSet.addOpt(new Opt(null, "out-stats", "Output stats to this file ", 1).setDefaultValue("/dev/stdout"));
//        optSet.addOpt(new Opt(null, "out-calls-AB", "Output calls to this file (in AB format)", 1));
//        optSet.addOpt(new Opt(null, "out-calls-IUPAC", "Output calls to this file (in IUPAC format) - due to the limitations "
//                + "of this format indels will be lost", 1));
//        optSet.addOpt(new Opt(null, "out-k-mer-coverage", "Output numbers ok k-mers covering each variant-base to this file ", 1));
//        optSet.addOpt(new Opt(null, "out-k-mer-frequencies", "Output median frequencies of k-mers covering each variant-base to this file", 1));
//        optSet.addOpt(new Opt(null, "out-calls-evidence", "Output base-call evidence to this file. This output is an amalgamation of the two previous ones", 1));
//        optSet.addOpt(new Opt(null, "out-calls-with-evidence", "Output base-call with evidence to this file. This output is an amalgamation of the two previous ones", 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
        optSet.addOpt(new Opt('D', "debug", "Print additional info for debugging purposes"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    private void buildFreqMerMap(OptSet optSet) {
        map = new HashMap<>();
        kmerFilters = new ArrayList<>();
        String fastaFileName = (String) optSet.getOpt("f").getValueOrDefault();
        int k = (int) optSet.getOpt("k").getValueOrDefault();
//        HashMap<String, Sequence> sequences = shared.FastaReader.hashMapOfSequencesFromFasta(fastaFileName, null);
        BufferedReader bufferdReader = null;
//        HashMap<String, ArrayList<KmerLink>> nonUniqueLinks = new HashMap<>();
        try {
            String inputLine;
            String id = "";
            StringBuilder seqBuilder = new StringBuilder();
            if (fastaFileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(fastaFileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(fastaFileName)), READER_BUFFER_SIZE);
            }
            while ((inputLine = bufferdReader.readLine()) != null) {
                String line = inputLine.trim();
                if (line.startsWith(">")) {
                    if (seqBuilder.length() > 0) {
                        KmerFilter kmerFilter = new KmerFilter(id, new Sequence(id, seqBuilder.toString()), k, TOOL_NAME);
                        if (kmerFilters.contains(kmerFilter)) {
                            continue;
                        }
                        kmerFilters.add(kmerFilter);
                        kmerizeAndAddToMap(kmerFilter, k, map);
                        seqBuilder = new StringBuilder();
                    }
                    id = line.substring(1); //get rid of ">"
                } else {
                    seqBuilder.append(line);
                }
//              
            }
            if (id.isEmpty()) {
                System.err.println("Error reading FASTA [" + fastaFileName + "]. No identifier lines? Terminating... ");
                System.exit(1);
            } else {
                KmerFilter kmerFilter = new KmerFilter(id, new Sequence(id, seqBuilder.toString()), k, TOOL_NAME);
                if (!kmerFilters.contains(kmerFilter)) {
                    kmerFilters.add(kmerFilter);
                    kmerizeAndAddToMap(kmerFilter, k, map);
                }

            }
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + fastaFileName, TOOL_NAME);
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bufferdReader != null) {
                    bufferdReader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        Reporter.report("[INFO]", intFormat(kmerFilters.size()) + " unique input sequences", TOOL_NAME);
        Reporter.report("[INFO]", intFormat(map.size()) + " k-mer-links in map", TOOL_NAME);
//        Reporter.report("[INFO]", intFormat(map.size()) + " k-mer-links in map, purging non-unique ones", TOOL_NAME);

    }
    //    private static void kmerizeAndAddToMap(CharSequence sequence, int k, int snpSite) {

    /**
     * 
     * @param kmerFilter
     * @param k
     * @param map
     * @return set of canonical k-mers for a given sequence
     */
    private void kmerizeAndAddToMap(KmerFilter kmerFilter, int k, HashMap<String, ArrayList<KmerLink>> map) { //,
//            HashMap<String, ArrayList<KmerLink>> nonUniqueLinks) {
        //TWO PARENT SEQUENCES FOR EACH SNP
        
        TreeSet<String> canonicalKmers = new TreeSet();
        String sequenceString = kmerFilter.getSequence().getSequenceString();
        int offset = 0; //If padding shifts snpSite, store that offset here
        int maxKmer = sequenceString.length() - k + 1;

//            if (DEBUG) {
//                for (int i = 0; i < sequence.length(); i++) {
//                    if (i == snpSite) {
//                        System.err.print("|");
//                    } else {
//                        System.err.print("_");
//                    }
//                }
//                System.err.println();
//                System.err.println(sequence);
//            }
        int startAt = 0; ///Math.max(0, snpSite - k + 1);
        for (int i = startAt; i < maxKmer; i++) {
            CharSequence kmer = sequenceString.subSequence(i, i + k);
            String canonical = SequenceOps.getCanonical(kmer.toString());

            int pos = i + offset; //position in the original/padded MSA sequence
//                KmerLink kmerLink = new KmerLink(snpFilter, (parent == 1), pos, !kmer.equals(canonical));
            KmerLink kmerLink = new KmerLink(kmerFilter, pos);
            ArrayList<KmerLink> kmerLinks = map.get(canonical);
            if (kmerLinks == null) {
                kmerLinks = new ArrayList();
            }
            kmerLinks.add(kmerLink);
            map.put(canonical, kmerLinks);
            canonicalKmers.add(canonical);
            if (DEBUG) {
                for (int j = 0; j < i; j++) {
                    System.err.print(" ");
                }
                System.err.println(kmer);
            }
        }
        if (DEBUG) {
            System.err.println();
        }
        kmerFilter.setMaxUniqMers(canonicalKmers.size());
    }

    private ArrayList<String> threadKmersThroughMap(OptSet optSet, ArrayList<String> kmersFileNames) {
        ArrayList<String> samples = new ArrayList<>();
        int IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();
        int IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();

        BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
//            boolean stranded = false;
        int ioThreads = 1 + 1;
        ArrayList<Future<?>> ioFutures = new ArrayList<>(ioThreads);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(ioThreads, ioThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, kmersFileNames, true, REPORTING_SHIFT, TOOL_NAME, IN_BUFFER_SIZE);
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
//        if (!inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.InFormat.KMERS)) {
//            Reporter.report("[FATAL]", "Only k-mer sets accepted as input. Guessed format: " + inputReaderProducer.getGuessedInputFormat(), getClass().getSimpleName());
//        } else {
        //Start KmergerConsumerProducer and OutputWriterConsumer threads
        int threads = 1;
        final ExecutorService execService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(threads * 5);
//        for (int i = 0; i < threads; i++) {
        futures.add(execService.submit(new CallerConsumer(inputQueue, TOOL_NAME, samples, map, kmerFilters)));
//        }
        execService.shutdown();
        ioExecutorService.shutdown();

        try {

            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            execService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
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

        return samples;
    }

    private void outputFasta(String outputFasta) {
        if (outputFasta != null) {
            try {
                BufferedWriter out;
                if (outputFasta.endsWith(".gz")) {
                    out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFasta)), "UTF-8"), WRITER_BUFFER_SIZE);
                } else {
                    out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFasta), "UTF-8"), WRITER_BUFFER_SIZE);
                }
                for (KmerFilter snpFilter : kmerFilters) {
                    StringBuilder sb = new StringBuilder(">");
                    sb.append(snpFilter.getId()).append("_").append(snpFilter.getSequence().getId());
                    sb.append(System.lineSeparator()).append(snpFilter.getSequence().getUnpaddedSequenceString());
                    sb.append(System.lineSeparator()).append(">");
                    out.write(sb.toString());
                    out.newLine();
                }
                out.flush();
            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }
            Reporter.report("[INFO]", "FASTA output sent to " + outputFasta, TOOL_NAME);
        }
    }

    private String intFormat(int value) {
        return NumberFormat.getIntegerInstance().format(value);
    }
    
    private void reportResults(ArrayList<String> samples, String outFile) {
        if (outFile != null) {
            int snpsGenotyped = 0;
            int indelsLost = 0;
            try {
                BufferedWriter out;
                if (outFile.endsWith(".gz")) {
                    out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile)), "UTF-8"), WRITER_BUFFER_SIZE);
                } else {
                    out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"), WRITER_BUFFER_SIZE);
                }

                StringBuilder header = new StringBuilder("SeqId");
                for (String sample : samples) {
                    header.append(DELIMITER).append(sample);
                }
                out.write(header.toString());
                out.newLine();

                
                for (KmerFilter kmerFilter : kmerFilters) {
                        snpsGenotyped++;
                        StringBuilder sb = new StringBuilder();
                        sb.append(kmerFilter.getId());
//                sb.append(snpFilter.getSequence1().getId()).append("__").append(snpFilter.getSequence2().getId());
//                        sb.append(DELIMITER).append(kmerFilter.getMers1()).append(DELIMITER);
//                        sb.append(DELIMITER).append(kmerFilter.getMersParent2());
                        for (String sample : samples) {
                            sb.append(DELIMITER).append(kmerFilter.getKmerFilterStats(sample).getMedianFrequency());                            
                            sb.append(DELIMITER).append(kmerFilter.getKmerFilterStats(sample).getCoverage());                            
                            sb.append(DELIMITER).append(kmerFilter.getKmerFilterStats(sample).getCoverageRatio());                            
                        }
//                    System.out.println(sb);
                        out.write(sb.toString());
                        out.newLine();
                }
                out.flush();
//            System.err.println(snpFilter.getSequence1().getId() + DELIMITER + snpFilter.getSequence2().getId()+" CALL: "+snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor));
//            System.err.println(snpFilter.getMedian1() + " <-- "+Arrays.toString(snpFilter.getMers1()));
//            System.err.println(snpFilter.getMedian2() + " <-- "+Arrays.toString(snpFilter.getMers2()));
            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }
//            Reporter.report("[INFO]", snpsGenotyped + " out of " + snpFilters.size() + " input SNVs genotyped and printed to " + outFile, TOOL_NAME);
        }
    }
}
