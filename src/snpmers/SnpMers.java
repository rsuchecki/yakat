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
package snpmers;

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
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListMap;
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
import shared.BaseCall;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;
import shared.Sequence;
import shared.StdRedirect;

/**
 *
 * @author rad
 */
public class SnpMers {

    private final static String DELIMITER = "\t";
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int REPORTING_SHIFT = 6; //input progress reporting every READER_BUFFER_SIZE*REPORTING_FACTOR records
    private final int READER_BUFFER_SIZE = 8192;
    private final int WRITER_BUFFER_SIZE = 8192;
    private String allele1 = "Al1";
    private String allele2 = "Al2";
//    private String parent1 = null;
//    private String parent2 = null;
    private ConcurrentSkipListMap<String, ArrayList<KmerLink>> map;
    private ArrayList<SnpFilter> snpFilters;

    private boolean DEBUG = false;

    private enum OutFmt {
        IUPAC, AB, SLASH, FREQ, COV, EVIDENCE, CALL_WITH_EVIDENCE, NEXUS;
    }

//    }
    public SnpMers(String[] args, String callerName, String toolName) {
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
        buildSnpMerMap(optSet);
//        removeSnpMersPoorlyCoveredInParents(optSet); //DONE AFTER FALSE POSITIVE FILTERING

        
        if (optSet.getOpt("F").isUsed()) {
            ArrayList<String> filteringKmersFileNames = (ArrayList<String>) optSet.getOpt("F").getValues();
            int minTotal = (int) optSet.getOpt("filter-min-k-mer-frequency-sum").getValueOrDefault();
            int minMinor = (int) optSet.getOpt("filter-min-k-mer-frequency-minor").getValueOrDefault();
            double minCoverage = (double) optSet.getOpt("filter-min-snp-coverage").getValueOrDefault();
            double maxError = (double) optSet.getOpt("filter-max-coverage-error").getValueOrDefault();
            ArrayList<String> filteringNames = threadKmersThroughMap(optSet, filteringKmersFileNames,
                minTotal, minMinor, minCoverage, maxError);
            removeFalsePositives(filteringNames);
        }
        int minTotal = (int) optSet.getOpt("min-k-mer-frequency-sum").getValueOrDefault();
        int minMinor = (int) optSet.getOpt("min-k-mer-frequency-minor").getValueOrDefault();
        double minCoverage = (double) optSet.getOpt("min-snp-coverage").getValueOrDefault();
        double maxError = (double) optSet.getOpt("max-coverage-error").getValueOrDefault();
        ArrayList<String> kmersFileNames = (ArrayList<String>) optSet.getOpt("K").getValues();
        ArrayList<String> sampleNames = threadKmersThroughMap(optSet, kmersFileNames, minTotal, minMinor, minCoverage, maxError);

        //Before reporting invalidate SNPs which do not meet reporting thresholds
        invalidateSnpsUnderThresholds(optSet);
        //Out
        outputFasta((String) optSet.getOpt("out-fasta").getValueOrDefault());
        outputNexus(sampleNames, (String) optSet.getOpt("out-distances-nexus").getValueOrDefault());
        reportResults(sampleNames, (String) optSet.getOpt("out-calls").getValueOrDefault(), optSet, OutFmt.SLASH);
        reportResults(sampleNames, (String) optSet.getOpt("out-calls-AB").getValueOrDefault(), optSet, OutFmt.AB);
        reportResults(sampleNames, (String) optSet.getOpt("out-calls-IUPAC").getValueOrDefault(), optSet, OutFmt.IUPAC);
        reportResults(sampleNames, (String) optSet.getOpt("out-k-mer-coverage").getValueOrDefault(), optSet, OutFmt.COV);
        reportResults(sampleNames, (String) optSet.getOpt("out-k-mer-frequencies").getValueOrDefault(), optSet, OutFmt.FREQ);
        reportResults(sampleNames, (String) optSet.getOpt("out-calls-evidence").getValueOrDefault(), optSet, OutFmt.EVIDENCE);
        reportResults(sampleNames, (String) optSet.getOpt("out-calls-with-evidence").getValueOrDefault(), optSet, OutFmt.CALL_WITH_EVIDENCE);
//        threadKmersThroughMap_BAK(optSet);
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("[A] Take a list of pre-defined SNPs and the associated sequences (e.g. LNISKS calls). "
            + "[B] Take a set of k-mers per sample of interest. "
            + "[C] Record frequencies of k-mers overlapping the input SNPs. "
            + "[D] For each pre-defined SNP, genotype each sample of interest");

        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt(null, "sample-ids", "Space separated sample identifiers which form the prefices of the input FASTA identifiers")
//            .setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE).setRequired(true));
        optSet.addOpt(new Opt('k', "k-mer-length", "It must match the size of the k-mers used to query the LNISKS snps", 1).setRequired(true).setMinValue(3).setMaxValue(255));
        optSet.addOpt(new Opt('s', "lnisks-snps", "File containing the table of SNPs called by LNISKS", 1).setRequired(true));
//        optSet.addOpt(new Opt('f', "lnisks-fasta", "The (msa) FASTA file matching the SNP information", 1));
        optSet.addOpt(new Opt('K', "per-sample-k-mers", "A set (or sets) of k-mers to be threaded through the map of k-mer-links to SNPs").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('F', "filtering-k-mers", "One or more sets of k-mers, each from one homozygous cultivar. "
            + "If both alleles of a putative SNP are present in such a set the SNP will be discarded as a likely false-positive").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 8192, 128, 65535));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up", 64, 1, 256));

//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Cluster processing settings]");
//        optSet.addOpt(new Opt(null, "min-samples-clustered", "Minimum number of samples in a cluster", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "min-seqs-clustered", "Minimum number of sequences in a cluster", 1).setMinValue(2).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "max-seqs-clustered", "Maximum number of sequences in a cluster", 1).setMinValue(2).setDefaultValue(1000));
        int footId = 1;
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Filtering-out input SNPs supported by few unique k-mers]");
//        optSet.addOpt(new Opt(null, "filter-min-uniqmers-ratio", "Remove parental SNPs if at least one allele is supported by fewer than <arg> ratio of unique"
//                + "k-mers to max possible k-mers", 1).setMinValue(0.01).setDefaultValue(0.75).addFootnote(footId, "There can be up to min(k,distance of snp from the end of the sequence) unique k-mers overlapping a SNP."));
//
//        footId++;
        String footText = "k-mer frequency corresponding to a SNP allele is obtained by taking a median of frequencies of all k-mers overlapping the underlying SNP";
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Filtering-out potential false-positive input SNPs]");
        optSet.addOpt(new Opt(null, "filter-min-k-mer-frequency-sum", "Minimum frequency of filtering k-mers "
            + "which overlap with a SNP (minor+major allele)", 1).setMinValue(1).setDefaultValue(2).addFootnote(footId, footText));
        optSet.addOpt(new Opt(null, "filter-min-k-mer-frequency-minor", "Minimum frequency of filtering k-mers "
            + "which support the minor allele", 1).setMinValue(1).setDefaultValue(1).addFootnote(footId, footText));
//        optSet.addOpt(new Opt(null, "filter-min-overlapping-k-mers", "At least <arg> filtering k-mers must overlap a locus (for each allele)", 1).setMinValue(1).setDefaultValue(1));
        optSet.addOpt(new Opt(null, "filter-min-snp-coverage", "At least <arg> fraction of unique snp-covering k-mers must be present in a filtering dataset (for each allele)", 1).setDefaultValue(0.75));
        optSet.addOpt(new Opt(null, "filter-max-coverage-error", "Maximum fraction of unique snp-covering k-mers that will be ignored as errorneous assignment and not considered for filtering", 1).setMinValue(0.00).setDefaultValue(0.05));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Variant calling and reporting]");
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-sum", "Minimum frequency of k-mers which overlap "
            + "with a SNP (minor+major allele)", 1).setMinValue(1).setDefaultValue(3).addFootnote(footId, footText));
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-minor", "Minimum frequency of k-mers "
            + "which support the minor allele", 1).setMinValue(1).setDefaultValue(3).addFootnote(footId, footText));
//        optSet.addOpt(new Opt(null, "min-overlapping-k-mers", "At least <arg> k-mers must overlap a locus (for each allele)", 1).setMinValue(1).setDefaultValue(1));
        optSet.addOpt(new Opt(null, "min-snp-coverage", "At least <arg> fraction of unique snp-covering k-mers must be present in a genotyped dataset (for each allele)", 1).setDefaultValue(0.9));
        optSet.addOpt(new Opt(null, "max-coverage-error", "Maximum fraction of unique snp-covering k-mers that will be ignored as errorneous assignment and not considered for genotyping", 1).setMinValue(0.00).setDefaultValue(0.05));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Reporting]");
        optSet.addOpt(new Opt(null, "min-samples-called", "Only report loci for which a base call was made for at least <arg> samples", 1).setMinValue(1).setDefaultValue(1));
        optSet.addOpt(new Opt(null, "min-samples-called-hom", "Only report loci for which a homozygous base  call was made for at least <arg> samples", 1).setMinValue(1).setDefaultValue(0));
        optSet.addOpt(new Opt(null, "min-samples-called-het", "Only report loci for which a heterozygous base call was made for at least <arg> samples", 1).setMinValue(0).setDefaultValue(0));
        optSet.addOpt(new Opt(null, "max-samples-called", "Only report loci for which a base call was made for at most <arg> samples", 1).setMinValue(1));
        optSet.addOpt(new Opt(null, "max-samples-called-hom", "Only report loci for which a homozygous base  call was made for at most <arg> samples", 1).setMinValue(0));
        optSet.addOpt(new Opt(null, "max-samples-called-het", "Only report loci for which a heterozygous base call was made for at most <arg> samples", 1).setMinValue(0));
//
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime and output settings]");
        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(1, "Multi-threading currently only implemented for initial snpMer map building"));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
        optSet.addOpt(new Opt(null, "out-fasta", "Output relevant sequences to this file (in FASTA format)", 1));
        optSet.addOpt(new Opt(null, "out-calls", "Output calls to this file ", 1).setDefaultValue("/dev/stdout"));
        optSet.addOpt(new Opt(null, "out-calls-AB", "Output calls to this file (in AB format)", 1));
        optSet.addOpt(new Opt(null, "out-calls-IUPAC", "Output calls to this file (in IUPAC format) - due to the limitations "
            + "of this format indels will be lost", 1));
        optSet.addOpt(new Opt(null, "out-k-mer-coverage", "Output numbers ok k-mers covering each variant-base to this file ", 1));
        optSet.addOpt(new Opt(null, "out-k-mer-frequencies", "Output median frequencies of k-mers covering each variant-base to this file", 1));
        optSet.addOpt(new Opt(null, "out-calls-evidence", "Output base-call evidence to this file. This output is an amalgamation of the two previous ones", 1));
        optSet.addOpt(new Opt(null, "out-calls-with-evidence", "Output base-call with evidence to this file. This output is an amalgamation of the two previous ones", 1));
        optSet.addOpt(new Opt(null, "out-distances-nexus", "Output distances between called samples in NEXUS format", 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
        optSet.addOpt(new Opt('D', "debug", "Print additional info for debugging purposes"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    private void buildSnpMerMap(OptSet optSet) {
        map = new ConcurrentSkipListMap<>();
        snpFilters = new ArrayList<>();

        String snpsFileName = (String) optSet.getOpt("s").getValueOrDefault();
        int k = (int) optSet.getOpt("k").getValueOrDefault();
        int threads = (int) optSet.getOpt("t").getValueOrDefault();

        int IN_Q_CAPACITY = threads; //(int) optSet.getOpt("Q").getValueOrDefault();
        int BUFFER_SIZE = 256; //(int) optSet.getOpt("U").getValueOrDefault();
        BlockingQueue snpFiltersQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
        final ExecutorService execService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
        for (int i = 0; i < threads; i++) {
            futures.add(execService.submit(new MapPopulatorConsumer(map, snpFiltersQueue, k, TOOL_NAME)));
        }
        ArrayList<SnpFilter> snpFiltersBuffer = new ArrayList<>(BUFFER_SIZE);

//        HashMap<String, Sequence> sequences = shared.FastaReader.hashMapOfSequencesFromFasta(fastaFileName, null);
        BufferedReader bufferdReader = null;
//        HashMap<String, ArrayList<KmerLink>> nonUniqueLinks = new HashMap<>();
        try {
            String inputLine;
            if (snpsFileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(snpsFileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(snpsFileName)), READER_BUFFER_SIZE);
            }
            Pattern spliPattern = Pattern.compile("\t");

//            SnpFilter previous = null;
//            String prevS1 = "";
//            String prevS2 = "";
            while ((inputLine = bufferdReader.readLine()) != null) {
                if (!inputLine.startsWith("ClusterId")) {
                    String line = inputLine.trim();
                    String[] toks = spliPattern.split(line);
//                    String[] toks = line.split(DELIMITER);
                    String clusterId = toks[0];
                    String id1 = toks[2];
                    String id2 = toks[4];
//                Sequence parent1 = sequences.get(clusterId + "_" + id1);
//                Sequence parent2 = sequences.get(clusterId + "_" + id2);
                    if (toks.length > 12) {
                        Reporter.report("[WARNING]", "Ignoring " + clusterId + " (" + toks.length + " cols)", TOOL_NAME);
                    } else {
//                        if (parent1 == null) {
//                            if (toks[2].contains("_")) {
//                                parent1 = toks[2].substring(0, toks[2].indexOf('_', 2));
//                                parent2 = toks[4].substring(0, toks[4].indexOf('_', 2));
//                            } else {
//                                parent1 = toks[2];
//                                parent2 = toks[4];
//                            }
//                        }
                        //If SNPs were called from multiple samples it is possible that 
                        //a SNP at a given position and involving the same nucleotides
                        //was called more than once. 
                        int snpPosition = Integer.parseInt(toks[6]) - 1;
//                        String s1 = toks[7].compareTo(toks[8]) <= 0 ? toks[7] : toks[8];
//                        String s2 = toks[7].compareTo(toks[8]) <= 0 ? toks[8] : toks[7];
//                        if (previous != null && previous.getClusterId().equals(clusterId) && previous.getSnpPosition0() == snpPosition
//                                && s1.equals(prevS1) && s2.equals(prevS2)) {
//                            continue;
//                        }
                        SnpFilter snpFilter = new SnpFilter(clusterId, new Sequence(id1, toks[10]), new Sequence(id2, toks[11]), snpPosition, TOOL_NAME);
                        if (snpFilters.contains(snpFilter)) {
                            continue;
                        }
                        snpFilters.add(snpFilter);
                        //                        kmerizeAndAddToMap(snpFilter, k, map, nonUniqueLinks);

                        snpFiltersBuffer.add(snpFilter);
                        if (snpFiltersBuffer.size() >= BUFFER_SIZE) {
                            snpFiltersQueue.put(snpFiltersBuffer);
                            snpFiltersBuffer = new ArrayList<>(BUFFER_SIZE);
                        }

//                        kmerizeAndAddToMap(snpFilter, k, map);
//                        previous = snpFilter;
//                        prevS1 = s1;
//                        prevS2 = s2;
                    }
                }
            }
            snpFiltersQueue.put(snpFiltersBuffer);
            snpFiltersQueue.put(new ArrayList());
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + snpsFileName, TOOL_NAME);
        } catch (IOException | InterruptedException ex) {
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

        execService.shutdown();

        try {

            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            execService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
//        for (Message fm : finalMessages) {
//            Reporter.report(fm.getLevel().toString(), fm.getBody(), fm.getCaller());
//        }

        Reporter.report("[INFO]", intFormat(snpFilters.size()) + " unique variant positions to be genotyped", TOOL_NAME);
        Reporter.report("[INFO]", intFormat(map.size()) + " k-mer-links in map", TOOL_NAME);
//        Reporter.report("[INFO]", intFormat(map.size()) + " k-mer-links in map, purging non-unique ones", TOOL_NAME);

//PURGE NON-UNIQUE k-mers (these can only be non-unique due to merging of non conflicting, clustered seeds from vsearch)
//        Iterator<Map.Entry<CharSequence, ArrayList<KmerLink>>> it = map.entrySet().iterator();
//        while (it.hasNext()) {
//            Map.Entry<CharSequence, ArrayList<KmerLink>> next = it.next();
//            ArrayList<KmerLink> kmerLinks = next.getValue();
//            if (kmerLinks.size() > 2) {
////                System.err.println(kmerLink.getParentSequence().getId()+" "+kmerLink.getParentSequence().getSequenceString()+" "+kmerLink.getStartPosition());
//                it.remove();
//            } else if (kmerLinks.size() == 2) {
//                KmerLink kLink0 = kmerLinks.get(0);
//                KmerLink kLink1 = kmerLinks.get(1);
//                String clusterId0 = kLink0.getSnpFilter().getClusterId();
//                String clusterId1 = kLink1.getSnpFilter().getClusterId();
//                int pos0 = 0;
//                int pos1 = 0;
//                if (clusterId0.contains(":")) { //non-LNISKS, mapping derived input
//                    String[] split = clusterId0.split(":");
//                    clusterId0 = split[0];
//                    pos0 = Integer.parseInt(split[1]);
//                }
//                if (clusterId1.contains(":")) {
//                    String[] split = clusterId1.split(":");
//                    clusterId1 = split[0];
//                    pos1 = Integer.parseInt(split[1]);
//                }
//                if (clusterId0.equals(clusterId1)) {
//                    int snpPosition0 = kLink0.getSnpFilter().getSnpPosition0();
//                    int snpPosition1 = kLink1.getSnpFilter().getSnpPosition0();
//                    int left = Math.min(snpPosition0, snpPosition1);
//                    int right = Math.max(snpPosition0, snpPosition1);
//
//                    if ((pos0 == pos1 && right - left + 1 <= k) || (pos0 != pos1 && Math.max(pos0, pos1) - Math.min(pos0, pos1) + 1 <= k)) {
//                        kLink0.getSnpFilter().incrementMerCount(kLink0.isParentOne());
//                        kLink1.getSnpFilter().incrementMerCount(kLink1.isParentOne());
//                    } else {
//                        it.remove();
//                    }
//                }
////            } else {
////                kmerLink.getSnpFilter();                
//            } else if (kmerLinks.size() < 2) {
//
//                //keep a tally of unique k-mers overlapping each SNP so that in genotyping this can be an upper bound required for a call 
//                for (KmerLink kmerLink : kmerLinks) {
//                    kmerLink.getSnpFilter().incrementMerCount(kmerLink.isParentOne());
//                }
//            }
//        }
        //INVESTIGATING CASES WHERE NON-UNIQUE k-mers HAVE BEEN OBSERVED 
//        Iterator<Map.Entry<String, ArrayList<KmerLink>>> iterator = nonUniqueLinks.entrySet().iterator();
//        while (iterator.hasNext()) {
//            Map.Entry<String, ArrayList<KmerLink>> next = iterator.next();
//            String key = next.getKey();
//            String[] toks = key.split(",");
//            if (!toks[0].equals(toks[1])) { //IGNORE MOST CASES FOR NOW
//                System.err.print(key + " ");
//                ArrayList<KmerLink> links = next.getValue();
//                Collections.sort(links);
//                int previous = -1;
//                for (KmerLink link : links) {
//                    int startPosition = link.getStartPosition();
//                    if (startPosition != previous) { //IGNORE MOST CASES FOR NOW
////                        System.err.print(startPosition + " ");
//                    }
//                    previous = startPosition;
//                }
//                System.err.println();
//                
//            }                       
//        }
//        Reporter.report("[INFO]", intFormat(map.size()) + " unique k-mer-links in map", TOOL_NAME);
    }
    //    private static void kmerizeAndAddToMap(CharSequence sequence, int k, int snpSite) {

    private ArrayList<String> threadKmersThroughMap(OptSet optSet, ArrayList<String> kmersFileNames,
        int minTotal, int minMinor, double minCoverage, double maxError) {
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
        int threads = 1;
        final ExecutorService execService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(threads * 5);
//        for (int i = 0; i < threads; i++) {
        futures.add(execService.submit(new CallerConsumer(inputQueue, TOOL_NAME, samples, map, snpFilters,
            minTotal, minMinor, minCoverage, maxError)));
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

    private void removeFalsePositives(ArrayList<String> filters) {
        Reporter.report("[INFO]", "Now threading false-positive filtering k-mer-s through k-mer-links-to-snps map...", TOOL_NAME);
        Iterator<SnpFilter> it = snpFilters.iterator();
        int invalid = 0;
        while (it.hasNext()) {
            SnpFilter snpFilter = it.next();
            for (String sample : filters) {
                //sb.append(DELIMITER).append(snpFilter.getBaseCall(sample));
                //IF A  HET WAS CALLED

                if (!snpFilter.getBaseCall(sample).getCallString().matches("[ACGTNacgtn-]")) {
//                    it.remove();                    
                    snpFilter.setInvalid();
                    invalid++;
                    if (DEBUG) {
                        Reporter.report("[INFO]", "SNV filtered-out " + snpFilter.getClusterId() + " "
                            + snpFilter.getSequence1().getId() + " " + snpFilter.getSequence2().getId() + " at " + snpFilter.getSnpPosition0(), TOOL_NAME);
//                    System.err.println("Call: " + snpFilter.getBaseCall(sample));
                    }
                    break; //ensures we only invalidate a snp once
                }
            }
        }
        Iterator<SnpFilter> it2 = snpFilters.iterator();
        SnpFilter previous = null;
        int invalidAdjacent = 0;
        while (it2.hasNext()) {
            SnpFilter snpFilter = it2.next();
            if (previous != null && snpFilter.getClusterId().equals(previous.getClusterId())) {
                if (!previous.isValid() && snpFilter.isValid() || previous.isValid() && !snpFilter.isValid()) {
//                    System.err.println("Only one of two SNVs filtered out");
                    previous.setInvalid();
                    snpFilter.setInvalid();
                    invalidAdjacent++;
                    invalid++;
                }
            }
            previous = snpFilter;
        }
        Reporter.report("[INFO]", "Filtered-out " + invalid + " likely fasle positives SNPs, including " + invalidAdjacent + " due to adjacent SNP being filtered out", TOOL_NAME);
    }

    private void removeSnpMersPoorlyCoveredInParents(OptSet optSet) {
        int k = (int) optSet.getOpt("k").getValueOrDefault();
        double minUniqRatio = (double) optSet.getOpt("filter-min-uniqmers-ratio").getValueOrDefault();
        Iterator<SnpFilter> iterator = snpFilters.iterator();
        int invalid = 0;
        System.err.println("P1\tP2\tID\tpos\tpos_rev");
        while (iterator.hasNext()) {
            SnpFilter snpFilter = iterator.next();
            int uniqeMersParent1 = snpFilter.getMersCountParent1();
//            if (snpFilter.getBase1() == '-') {
//                uniqeMersParent1++;
//            }
            int uniqeMersParent2 = snpFilter.getMersCountParent2();
//            if (snpFilter.getBase2() == '-') {6
//                uniqeMersParent2++;
//            }
//            System.err.println();
            int lengthUnpadded1 = snpFilter.getSequence1().getLengthUnpadded();
            int lengthUnpadded2 = snpFilter.getSequence2().getLengthUnpadded();
            int positionUnpadded1 = snpFilter.getSnpPosition0UnpaddedSeq1();
            int positionUnpadded2 = snpFilter.getSnpPosition0UnpaddedSeq2();

//            int maxUniq1 = Math.min(positionUnpadded1 + 1, lengthUnpadded1 - positionUnpadded1);
//            int maxUniq2 = Math.min(positionUnpadded2 + 1, lengthUnpadded2 - positionUnpadded2);
            int leftDelatK1 = k - positionUnpadded1 - 1;
            int rightDelatK1 = k - (lengthUnpadded1 - positionUnpadded1);
            int maxUniq1 = k - Math.max(leftDelatK1, 0) - Math.max(rightDelatK1, 0);
            int leftDelatK2 = k - positionUnpadded2 - 1;
            int rightDelatK2 = k - (lengthUnpadded2 - positionUnpadded2);
            int maxUniq2 = k - Math.max(leftDelatK2, 0) - Math.max(rightDelatK2, 0);

            maxUniq1 = k < maxUniq1 ? k : maxUniq1;
            maxUniq2 = k < maxUniq2 ? k : maxUniq2;
            double uniqCov1 = (double) uniqeMersParent1 / maxUniq1;
            double uniqCov2 = (double) uniqeMersParent2 / maxUniq2;

            if (uniqCov1 < minUniqRatio || uniqCov2 < minUniqRatio) {
                System.err.println("Removing\t" + uniqeMersParent1 + " of " + maxUniq1 + "\t" + uniqeMersParent2 + " of " + maxUniq2 + "\t"
                    + snpFilter.getClusterId() + "\t" + (snpFilter.getSnpPosition0UnpaddedSeq1() + 1) + "\t"
                    + (snpFilter.getSequence1().getLength() - snpFilter.getSnpPosition0UnpaddedSeq2()) + "\t"
                    + snpFilter.getSequence1().getSequenceString() + "\t" + snpFilter.getSequence2().getSequenceString());
                snpFilter.setInvalid();
                invalid++;
            } else {
//                System.err.println("KEEP\t" + uniqeMersParent1 + " of " + maxUniq1 + "\t" + uniqeMersParent2 + " of " + maxUniq2+ "\t" 
//                        + snpFilter.getClusterId() + "\t" + (snpFilter.getSnpPosition0UnpaddedSeq1()+ 1) + "\t" 
//                        + (snpFilter.getSequence1().getLength() - snpFilter.getSnpPosition0UnpaddedSeq2())+ "\t"
//                        + snpFilter.getSequence1().getSequenceString() + "\t" +snpFilter.getSequence2().getSequenceString());
//                System.out.print("");
            }
        }
        Reporter.report("[INFO]", "Filtered-out " + invalid + " SNPs due to unique k-mer coverage under " + minUniqRatio, TOOL_NAME);
    }

    private void invalidateSnpsUnderThresholds(OptSet optSet) {
        int minCalled = (int) optSet.getOpt("min-samples-called").getValueOrDefault();
        int minCalledHom = (int) optSet.getOpt("min-samples-called-hom").getValueOrDefault();
        int minCalledHet = (int) optSet.getOpt("min-samples-called-het").getValueOrDefault();

        Integer maxCalled = (Integer) optSet.getOpt("max-samples-called").getValueOrDefault(null);
        Integer maxCalledHom = (Integer) optSet.getOpt("max-samples-called-hom").getValueOrDefault(null);
        Integer maxCalledHet = (Integer) optSet.getOpt("max-samples-called-het").getValueOrDefault(null);
        
        for (SnpFilter snpFilter : snpFilters) {
            int callsHom = snpFilter.getHomozygousCalls();
            int callsHet = snpFilter.getHeterozygousCalls();
            int calls = callsHet + callsHom;
            if(calls < minCalled || callsHom < minCalledHom || callsHet < minCalledHet) {
                snpFilter.setInvalid();
                continue;
            }
            if(maxCalled != null && calls > maxCalled) {
                snpFilter.setInvalid();
                continue;                
            }
            if(maxCalledHom != null && callsHom > maxCalledHom) {
                snpFilter.setInvalid();
                continue;                
            }
            if(maxCalledHet != null && callsHet > maxCalledHet) {
                snpFilter.setInvalid();
                continue;                
            }
            Collection<BaseCall> values = snpFilter.getSnpCalls().values();
            int aa = 0;
            int bb = 0;
            
            //TODO make this optional
            for (BaseCall value : values) {
                String callAB = value.getCallAB(snpFilter.getBase1(), snpFilter.getBase2());
                if("AA".equals(callAB)) {
                  aa++;
                } else if ("BB".equals(callAB)) {
                    bb++;
                }
            }
            if(aa==0 || bb==0) {
                snpFilter.setInvalid();
            }
        }

    }

    private void reportResults(ArrayList<String> samples, String outFile, OptSet optSet, OutFmt fmt) {

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

                StringBuilder header = new StringBuilder("Ref:Pos");
                if (fmt == OutFmt.AB) {
                    header.append(DELIMITER).append("Phenotype");
                }
                header.append(DELIMITER).append(allele1);
                if (fmt != OutFmt.AB) {
                    header.append(DELIMITER).append("Cov1");
                }
                header.append(DELIMITER).append(allele2);
                if (fmt != OutFmt.AB) {
                    header.append(DELIMITER).append("Cov2");
                }
                for (String sample : samples) {
                    header.append(DELIMITER).append(sample).append("_call");
                    if (fmt == OutFmt.CALL_WITH_EVIDENCE) {
                        header.append(DELIMITER).append(sample).append("_cov1");
                        header.append(DELIMITER).append(sample).append("_cov2");
                        header.append(DELIMITER).append(sample).append("_freq1");
                        header.append(DELIMITER).append(sample).append("_freq2");
                    }
                }
//            System.out.println(header);
                out.write(header.toString());
                out.newLine();

                for (SnpFilter snpFilter : snpFilters) {
                    if (snpFilter.isValid()) {
                        snpsGenotyped++;
                        StringBuilder sb = new StringBuilder();
                        sb.append(snpFilter.getClusterId());
//                sb.append(snpFilter.getSequence1().getId()).append("__").append(snpFilter.getSequence2().getId());
                        sb.append(":").append(snpFilter.getSnpPosition0() + 1).append(DELIMITER);

                        if (fmt == OutFmt.AB) {
                            sb.append("All").append(DELIMITER);
                            sb.append("AA").append(DELIMITER);//.append(snpFilter.getMersParent1()).append(DELIMITER);
                            sb.append("BB");//.append(DELIMITER).append(snpFilter.getMersParent1());
                        } else {
                            if (fmt == OutFmt.IUPAC && snpFilter.isIndel()) {
                                --snpsGenotyped;
                                indelsLost++;
                                continue;
                            }
                            sb.append(snpFilter.getBase1()).append(DELIMITER).append(snpFilter.getMersCountParent1()).append(DELIMITER);
                            sb.append(snpFilter.getBase2()).append(DELIMITER).append(snpFilter.getMersCountParent2());
                        }
                        for (String sample : samples) {
                            sb.append(DELIMITER);
                            switch (fmt) {
                                case AB:
                                    sb.append(snpFilter.getBaseCall(sample).getCallAB(snpFilter.getBase1(), snpFilter.getBase2()));
                                    break;
                                case IUPAC:
                                    sb.append(snpFilter.getBaseCall(sample).getCallIUPAC());
                                    break;
                                case SLASH:
                                    sb.append(snpFilter.getBaseCall(sample).getCallString());
                                    break;
                                case COV:
                                    sb.append(snpFilter.getBaseCall(sample).getCoverageString());
                                    break;
                                case FREQ:
                                    sb.append(snpFilter.getBaseCall(sample).getFreqString());
                                    break;
                                case EVIDENCE:
                                    sb.append(snpFilter.getBaseCall(sample).getEvidence());
                                    break;
                                case CALL_WITH_EVIDENCE:
                                    sb.append(snpFilter.getBaseCall(sample).getCallString());
                                    sb.append(DELIMITER);
                                    sb.append(snpFilter.getBaseCall(sample).getCov1());
                                    sb.append(DELIMITER);
                                    sb.append(snpFilter.getBaseCall(sample).getCov2());
//                                    sb.append(snpFilter.getBaseCall(sample).getCoverageString());
                                    sb.append(DELIMITER);
                                    sb.append(snpFilter.getBaseCall(sample).getFreq1String());
                                    sb.append(DELIMITER);
                                    sb.append(snpFilter.getBaseCall(sample).getFreq2String());
//                                    sb.append(snpFilter.getBaseCall(sample).getFreqString());
                                    break;
                                default:
                                    break;
                            }
                        }
//                    System.out.println(sb);
                        out.write(sb.toString());
                        out.newLine();
                    }
                }
                out.close();
//            System.err.println(snpFilter.getSequence1().getId() + DELIMITER + snpFilter.getSequence2().getId()+" CALL: "+snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor));
//            System.err.println(snpFilter.getMedian1() + " <-- "+Arrays.toString(snpFilter.getMers1()));
//            System.err.println(snpFilter.getMedian2() + " <-- "+Arrays.toString(snpFilter.getMers2()));
            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }
            if (indelsLost != 0) {
                Reporter.report("[INFO]", indelsLost + " indels lost in the IUPAC output: " + outFile, TOOL_NAME);
            }
            Reporter.report("[INFO]", snpsGenotyped + " out of " + snpFilters.size() + " input SNVs genotyped and printed to " + outFile, TOOL_NAME);
        }
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
                for (SnpFilter snpFilter : snpFilters) {
                    if (snpFilter.isValid()) {
                        StringBuilder sb = new StringBuilder(">");
                        sb.append(snpFilter.getClusterId()).append("_").append(snpFilter.getSequence1().getId());
                        sb.append(":").append(snpFilter.getSnpPosition0UnpaddedSeq1() + 1);
                        sb.append(System.lineSeparator()).append(snpFilter.getSequence1().getUnpaddedSequenceString());
                        sb.append(System.lineSeparator()).append(">");
                        sb.append(snpFilter.getClusterId()).append("_").append(snpFilter.getSequence2().getId());
                        sb.append(":").append(snpFilter.getSnpPosition0UnpaddedSeq2() + 1);
                        sb.append(System.lineSeparator()).append(snpFilter.getSequence2().getUnpaddedSequenceString());
                        out.write(sb.toString());
                        out.newLine();
                    }
                }
                out.close();
            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }
            Reporter.report("[INFO]", "FASTA output sent to " + outputFasta, TOOL_NAME);
        }
    }

    private void outputNexus(ArrayList<String> samples, String outputNexus) {
        if (outputNexus != null) {
            try {
                BufferedWriter out;
                if (outputNexus.endsWith(".gz")) {
                    out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputNexus)), "UTF-8"), WRITER_BUFFER_SIZE);
                } else {
                    out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputNexus), "UTF-8"), WRITER_BUFFER_SIZE);
                }
                StringBuilder sb = new StringBuilder("#nexus");
                sb.append(System.lineSeparator()).append(System.lineSeparator());;
                sb.append("BEGIN Taxa;").append(System.lineSeparator());
                sb.append("DIMENSIONS ntax=").append(samples.size()).append(";").append(System.lineSeparator());
                sb.append("TAXLABELS").append(System.lineSeparator());
                for (int i = 0; i < samples.size(); i++) {
                    sb.append("[").append(i).append("] '").append(samples.get(i)).append("'").append(System.lineSeparator());
                }
                sb.append(";").append(System.lineSeparator());
                sb.append("END; [Taxa]").append(System.lineSeparator());
                sb.append(System.lineSeparator());
                sb.append("BEGIN Distances;").append(System.lineSeparator());
                sb.append("DIMENSIONS ntax=").append(samples.size()).append(";").append(System.lineSeparator());
                sb.append("FORMAT labels=left diagonal triangle=both;").append(System.lineSeparator());
                sb.append("MATRIX").append(System.lineSeparator());
                out.write(sb.toString());
                int[][] matches = new int[samples.size()][samples.size()];
                int[][] called = new int[samples.size()][samples.size()];

                for (SnpFilter snpFilter : snpFilters) {
                    if (snpFilter.isValid()) {
                        HashMap<String, BaseCall> snpCalls = snpFilter.getSnpCalls();
                        for (int i = 0; i < samples.size() - 1; i++) {
                            BaseCall call = snpCalls.get(samples.get(i));
                            if (call.isCalled()) {
                                for (int j = i + 1; j < samples.size(); j++) {
                                    BaseCall anotherCall = snpCalls.get(samples.get(j));
                                    if (anotherCall.isCalled()) {
                                        called[i][j]++;
                                        called[j][i]++;
                                        if (call.getCallString().equals(anotherCall.getCallString())) {
                                            matches[i][j]++;
                                            matches[j][i]++;
                                        }
                                    }
                                }
                            }

                        }

                    }
                }
                //DISTANCES COLLECTED, PRINT
                for (int i = 0; i < samples.size(); i++) {
                    sb = new StringBuilder();
                    sb.append("[").append(i).append("] '").append(samples.get(i)).append("'");
                    for (int j = 0; j < samples.size(); j++) {
                        if (i == j) {
                            sb.append(" ").append(0);
                        } else {
                            double distance = 1 - (double) matches[i][j] / (double) called[i][j];
                            sb.append(" ").append(distance);
                        }
                    }
                    out.write(sb.toString());
                    out.newLine();
                }

//                for (SnpFilter snpFilter : snpFilters) {                    
//                    if (snpFilter.isValid()) {
//                        
//                    }
//                }
                sb = new StringBuilder();
                sb.append(";").append(System.lineSeparator());
                sb.append("END; [Distances]").append(System.lineSeparator());
                out.write(sb.toString());
                out.close();

            } catch (UnsupportedEncodingException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            } catch (IOException e) {
                Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
            }

            Reporter.report(
                "[INFO]", "NEXUS distance matrix output sent to " + outputNexus, TOOL_NAME);
        }
    }

    private String intFormat(int value) {
        return NumberFormat.getIntegerInstance().format(value);
    }
}
