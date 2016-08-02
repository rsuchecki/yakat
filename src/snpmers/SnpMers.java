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
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import shared.FastaReader;
import shared.Reporter;
import shared.Sequence;
import shared.SequenceOps;
import shared.StdRedirect;

/**
 *
 * @author rad
 */
public class SnpMers {

    private final static String DELIMITER = "\t";
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int READER_BUFFER_SIZE = 8192;
    private final int WRITER_BUFFER_SIZE = 8192;
    private String parent1 = null;
    private String parent2 = null;
    private HashMap<CharSequence, KmerLink> map;
    ArrayList<SnpFilter> snpFilters;

    private boolean DEBUG = false;

//    }
    public SnpMers(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }
        if (optSet.getOpt("D").isUsed()) {
            DEBUG = true;
        }
        new StdRedirect(optSet, TOOL_NAME);

        buildSnpMerMap(optSet);
        filterOutFalsePositiveSnps(optSet);
        threadKmersThroughMap(optSet);

//        readAndProcessMSASequencesFromFasta(fileName, optSet);
//        readMSASequencesFromFasta(optSet.getPositionalOptsList().get(0).getValue());
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Given a list of NIKS-derived inter-parent SNPs and the corresponding "
            + "MSA FASTA, read k-mers from offspring samples and record frequencies of k-mers overlapping with "
            + "the input SNPs. These can then be used to call offsping base for a given parental SNP.");

        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt(null, "sample-ids", "Space separated sample identifiers which form the prefices of the input FASTA identifiers")
//            .setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE).setRequired(true));
        optSet.addOpt(new Opt('k', "k-mer-length", "It must match the size of the k-mers used to query the niks-snps", 1).setRequired(true).setMinValue(3).setMaxValue(255));
        optSet.addOpt(new Opt('s', "niks-snps", "File containing the table of SNPs called by NIKS", 1));
//        optSet.addOpt(new Opt('f', "niks-fasta", "The (msa) FASTA file matching the SNP information", 1));
        optSet.addOpt(new Opt('F', "filtering-k-mers", "One or more sets of k-mers, each from one homozygous cultivar. "
            + "If both alleles of a putative SNP are present in such a set the SNP will be discarded as a likely false-positive").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));

//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Cluster processing settings]");
//        optSet.addOpt(new Opt(null, "min-samples-clustered", "Minimum number of samples in a cluster", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "min-seqs-clustered", "Minimum number of sequences in a cluster", 1).setMinValue(2).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "max-seqs-clustered", "Maximum number of sequences in a cluster", 1).setMinValue(2).setDefaultValue(1000));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Variant calling and reporting]");
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-sum", "Minimum frequency of k-mers which overlap with a SNP (minor+major allele)", 1).setMinValue(1).setDefaultValue(5));
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-minor", "Minimum frequency of k-mers which support the minor allele", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "min-overlapping-k-mers", "At least <arg> k-mers must overlap a locus (for each allele)", 1).setMinValue(1).setDefaultValue(1));
        
////        optSet.addOpt(new Opt(null, "min-sequences-per-cluster", "Minimum number of sequences required for a cluster to be considered ", 1).setMinValue(2).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "min-inter-identity", "Minimum inter sample identity ", 1).setMinValue(0.0).setDefaultValue(0.95));
//        optSet.addOpt(new Opt(null, "max-inter-snps", "SNPs will be reported if at most <arg> inter-sample SNPs are called in a cluster", 1).setMinValue(0).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "max-intra-snps", "SNPs will be reported if at most <arg> intra-sample SNPs are called in a cluster", 1).setMinValue(0).setDefaultValue(1));
//        optSet.addOpt(new Opt(null, "max-indel-length", "Maximum length of an indel ", 1).setMinValue(0).setDefaultValue(1));
//        optSet.addOpt(new Opt(null, "min-indel-distance", "Treat '-' positions less than <arg> bases from sequence ends as padding not indels", 1).setMinValue(1).setDefaultValue(3));
//        optSet.addOpt(new Opt(null, "supress-inter-snps", "Do not report inter-sample SNPs"));
//        optSet.addOpt(new Opt(null, "supress-intra-snps", "Do not report intra-sample SNPs"));
//        optSet.addOpt(new Opt(null, "reverse-lex-order", "Variants will be reported in reverse lexicographical order of the sample/bulk identifiers"));
//
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Runtime and output settings]");
////        String threadsOrderNote = "Note that in multi-threaded mode the output lines order need not reflect the input order";
////        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(1, threadsOrderNote));
////        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 1024, 128, 32768));
////        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",64, 1, 256));
//////        optSet.addOpt(new Opt('u', "out-buffer-size", "Size of buffers put on out-queue ", 1024, 128, 32768));
//////        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writing-out",64, 1, 256));
//        optSet.addOpt(new Opt(null, "out-clusters-msa", "Output clustered sequences (for which SNPs were called) to <arg> MSA/FASTA file", 1));
//        optSet.addOpt(new Opt(null, "out-unclustered-fasta", "Output unclustered sequences to <arg> FASTA file", 1));
//        optSet.addOpt(new Opt(null, "out-unclustered-min-len", "Minimum length required to output an unclustered sequence", 1).setMinValue(1).setDefaultValue(100));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
////        String headerNote = "Can be useful for external parallization (print header once)";
////        optSet.addOpt(new Opt('H', "header-only", "Print header and exit").addFootnote(1, TOOL_NAME));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[A little bit of help]");
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
//        optSet.addOpt(new Opt('I', "iupac-codes-table", "Print the table of IUPAC nucleotide codes and exit"));
        optSet.addOpt(new Opt('D', "debug", "Print additional info for debugging purposes"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    private void buildSnpMerMap(OptSet optSet) {
        map = new HashMap<>();
        snpFilters = new ArrayList<>();

        String snpsFileName = (String) optSet.getOpt("s").getValueOrDefault();
        int k = (int) optSet.getOpt("k").getValueOrDefault();
//        HashMap<String, Sequence> sequences = shared.FastaReader.hashMapOfSequencesFromFasta(fastaFileName, null);
        BufferedReader bufferdReader = null;
        HashMap<String, ArrayList<KmerLink>> nonUniqueLinks = new HashMap<>();
        try {
            String inputLine;
            if (snpsFileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(snpsFileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(snpsFileName)), READER_BUFFER_SIZE);
            }
            while ((inputLine = bufferdReader.readLine()) != null) {
                if (!inputLine.startsWith("ClusterId")) {
                    String line = inputLine.trim();
                    String[] toks = line.split(DELIMITER);
                    String clusterId = toks[0];
                    String id1 = toks[2];
                    String id2 = toks[4];
//                Sequence parent1 = sequences.get(clusterId + "_" + id1);
//                Sequence parent2 = sequences.get(clusterId + "_" + id2);
                    if (toks.length > 12) {
                        Reporter.report("[WARNING]", "Ignoring " + clusterId + " (" + toks.length + " cols)", TOOL_NAME);
                    } else {
                        if(parent1 == null) {
                            parent1 = toks[2].substring(0, toks[2].indexOf('_'));
                            parent2 = toks[4].substring(0, toks[4].indexOf('_'));
                        }
                        int snpPosition = Integer.parseInt(toks[6]) - 1;
                        SnpFilter snpFilter = new SnpFilter(clusterId, new Sequence(id1, toks[10]), new Sequence(id2, toks[11]), snpPosition);
                        snpFilters.add(snpFilter);
                        kmerizeAndAddToMap(snpFilter, k, map, nonUniqueLinks);
                    }
                }
            }
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + snpsFileName, TOOL_NAME);
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
        Reporter.report("[INFO]", intFormat(map.size()) + " k-mer-links in map, purging non-unique ones", TOOL_NAME);
        //PURGE NON-UNIQUE k-mers (these can only be non-unique due to merging of non conflicting, clustered seeds from vsearch)
        Iterator<Map.Entry<CharSequence, KmerLink>> it = map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<CharSequence, KmerLink> next = it.next();
            KmerLink kmerLink = next.getValue();
            if (!kmerLink.isUnique()) {
//                System.err.println(kmerLink.getParentSequence().getId()+" "+kmerLink.getParentSequence().getSequenceString()+" "+kmerLink.getStartPosition());
                it.remove();
//            } else {
//                kmerLink.getSnpFilter();                
            }
        }
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
        Reporter.report("[INFO]", intFormat(map.size()) + " unique k-mer-links in map", TOOL_NAME);
//        it = map.entrySet().iterator();
//        while (iterator.hasNext()) {
//            Map.Entry<String, ArrayList<KmerLink>> next = iterator.next();
//            
//        }

    }

//    private static void kmerizeAndAddToMap(CharSequence sequence, int k, int snpSite) {
    private void kmerizeAndAddToMap(SnpFilter snpFilter, int k, HashMap<CharSequence, KmerLink> map,
        HashMap<String, ArrayList<KmerLink>> nonUniqueLinks) {
        //TWO PARENT SEQUENCES FOR EACH SNP
        for (int parent = 1; parent < 3; parent++) {
            Sequence parentSequence = snpFilter.getParentSequence(parent);
            String sequence = parentSequence.getSequenceString();
            int snpSite = snpFilter.getSnpPosition();
            int offset = 0; //If padding shifts snpSite, store that offset here
            int maxKmer = Math.min(sequence.length() - k + 1, snpSite + 1);

            if (DEBUG) {
                for (int i = 0; i < sequence.length(); i++) {
                    if (i == snpSite) {
                        System.err.print("|");
                    } else {
                        System.err.print("_");
                    }
                }
                System.err.println();
                System.err.println(sequence);
            }

            //CANNOT K-MERIZE GAPS INSERTED BY MSA, SO NEED TO REMOVE THEM AND ADJUST THE SNP POSITION ACCORDINGLY
            if (sequence.contains("-")) {
                StringBuilder unpaddedSeq = new StringBuilder();
                int unpaddedSnpPosition = snpSite;
                for (int i = 0; i < sequence.length(); i++) {
                    if (sequence.charAt(i) == '-') {
                        if (i < snpSite) {
                            --unpaddedSnpPosition;
                            offset++;
                        }
                    } else {
                        unpaddedSeq.append(sequence.charAt(i));
                    }
                }
                sequence = unpaddedSeq.toString();
                snpSite = unpaddedSnpPosition;
                maxKmer = Math.min(sequence.length() - k + 1, snpSite + 1);

                if (DEBUG) {
                    for (int i = 0; i < sequence.length(); i++) {
                        if (i == snpSite) {
                            System.err.print("|");
                        } else {
                            System.err.print("_");
                        }
                    }
                    System.err.println();
                    System.err.println(sequence);
                }
            }

            int startAt = Math.max(0, snpSite - k + 1);
            for (int i = startAt; i < maxKmer; i++) {
                CharSequence kmer = sequence.subSequence(i, i + k);
                CharSequence canonical = SequenceOps.getCanonical(kmer);

                int pos = i + offset; //position in the original/padded MSA sequence
                KmerLink kmerLink = new KmerLink(snpFilter, (parent == 1), pos, !kmer.equals(canonical));
                if (map.containsKey(canonical)) {
                    kmerLink.setUnique(false);
                }
                KmerLink put = map.put(canonical, kmerLink);
                if (put != null) {
                    String key = parentSequence.getId() + "," + put.getParentSequence().getId();

                    ArrayList<KmerLink> nonUniqList = nonUniqueLinks.getOrDefault(key, new ArrayList<KmerLink>());
                    if (!nonUniqList.contains(kmerLink)) {
                        nonUniqList.add(kmerLink);
                    }
                    if (!nonUniqList.contains(put)) {
                        nonUniqList.add(put);
                    }
                    nonUniqueLinks.put(key, nonUniqList);
//                    Reporter.report("[WARNING]", "Non-unique ["+parentSequence.getId()+","+put.getParentSequence().getId()+"] k-mer overlapping a SNP will be ignored: " + canonical, TOOL_NAME);
//                        Reporter.report("[ERROR]", "Unexpected issue of a k-mer already present in the map: " + canonical, TOOL_NAME);
                    if (DEBUG) {
                        System.err.println("kmer, canonical:");
                        System.err.println(kmer + ", " + canonical);
//                        System.err.println("Current:\n" + parentSequence.getId() + " at0=" + startPosition);
                        System.err.println("Current:\n" + parentSequence.getId() + " at0=" + pos);
                        System.err.println(parentSequence.getSequenceString());
                        System.err.println("Previous:\n" + put.getParentSequence().getId() + " at0=" + put.getStartPosition());
                        System.err.println(put.getParentSequence().getSequenceString());
                    }
                }

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

        }
    }

    private void filterOutFalsePositiveSnps(OptSet optSet) {
        ArrayList<String> kmersFileNames = (ArrayList<String>) optSet.getOpt("F").getValues();        
        Reporter.report("[INFO]", "Now threading false-positive filtering k-mer-s through k-mer-links-to-snps map...", TOOL_NAME);
//        ArrayList<String> samples = new ArrayList<>();
//        if (kmersFileNames.isEmpty()) {
//            kmersFileNames.add(null);
//        }
//
//        for (String kmersFileName : kmersFileNames) {
//            String sampleName = null;
//            BufferedReader bufferdReader = null;
//            try {
//                String inputLine;
//                if (kmersFileName == null) {
////                Reporter.report("[INFO]", "Input file(s) not specified, reading from stdin ", TOOL_NAME);
//                    bufferdReader = new BufferedReader(new InputStreamReader(System.in), READER_BUFFER_SIZE);
//                    if ((inputLine = bufferdReader.readLine()) != null) {
//                        sampleName = inputLine; //if stding first line is expected to be sample name
//                    }
//                } else if (kmersFileName.endsWith(".gz")) {
//                    InputStream gzipStream = new GZIPInputStream(new FileInputStream(kmersFileName), READER_BUFFER_SIZE);
//                    bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
//                    sampleName = kmersFileName;
//                } else {
//                    bufferdReader = new BufferedReader(new FileReader(new File(kmersFileName)), READER_BUFFER_SIZE);
//                    sampleName = kmersFileName;
//                }
//                Reporter.report("[INFO]", "Current sample: " + sampleName, TOOL_NAME);
//                while ((inputLine = bufferdReader.readLine()) != null) {
//                    String line = inputLine.trim();
//                    String[] toks = line.split(DELIMITER);
//                    if (toks.length == 1) {
//                        for (SnpFilter snpFilter : snpFilters) {
//                            snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor, minKmers, TOOL_NAME);
//                        }
//                        samples.add(sampleName);
//                        sampleName = line;
//                        Reporter.report("[INFO]", "Current sample: " + sampleName, TOOL_NAME);
//                    } else if (toks.length != 2) {
//                        Reporter.report("[ERROR]", "Unexpected input line: " + line, TOOL_NAME);
//                        System.exit(1);
//                    }
//                    KmerLink kmerLink = map.get(toks[0]);
//                    if (kmerLink != null) {
//                        boolean setMer = kmerLink.setMer(Short.parseShort(toks[1]));
//                        if (!setMer) {
//                            System.err.println("mer not set");
//                        }
////                    SnpFilter snpFilter = kmerLink.getSnpFilter();
////                    System.err.println(kmerLink.getParentSequence().getId()+"\t"+snpFilter.getSnpPosition()+"\t"+toks[1]);
//                    }
//                }
//                for (SnpFilter snpFilter : snpFilters) {
//                    snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor, minKmers, TOOL_NAME);
//                }
//                samples.add(sampleName);
//            } catch (FileNotFoundException ex) {
//                Reporter.report("[ERROR]", "File not found: " + kmersFileName, TOOL_NAME);
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            } finally {
//                try {
//                    if (bufferdReader != null) {
//                        bufferdReader.close();
//                    }
//                } catch (IOException ex) {
//                    ex.printStackTrace();
//                }
//            }
//        }
        Reporter.report("[INFO]", "Finished filtering -out likely fasle positives SNPs", TOOL_NAME);
    }
    
    private void threadKmersThroughMap(OptSet optSet) {
        ArrayList<String> kmersFileNames = (ArrayList<String>) optSet.getOpt("K").getValues();
//        String kmersFileName = (String) optSet.getOpt("K").getValueOrDefault();
        int minTotal = 5;
        int minMinor = 2;
        int minKmers = 6;
        Reporter.report("[INFO]", "Now threading input k-mer-s through k-mer-links-to-snps map...", TOOL_NAME);
        ArrayList<String> samples = new ArrayList<>();
        if (kmersFileNames.isEmpty()) {
            kmersFileNames.add(null);
        }

        for (String kmersFileName : kmersFileNames) {
            String sampleName = null;
            BufferedReader bufferdReader = null;
            try {
                String inputLine;
                if (kmersFileName == null) {
//                Reporter.report("[INFO]", "Input file(s) not specified, reading from stdin ", TOOL_NAME);
                    bufferdReader = new BufferedReader(new InputStreamReader(System.in), READER_BUFFER_SIZE);
                    if ((inputLine = bufferdReader.readLine()) != null) {
                        sampleName = inputLine; //if stding first line is expected to be sample name
                    }
                } else if (kmersFileName.endsWith(".gz")) {
                    InputStream gzipStream = new GZIPInputStream(new FileInputStream(kmersFileName), READER_BUFFER_SIZE);
                    bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
                    sampleName = kmersFileName;
                } else {
                    bufferdReader = new BufferedReader(new FileReader(new File(kmersFileName)), READER_BUFFER_SIZE);
                    sampleName = kmersFileName;
                }
                Reporter.report("[INFO]", "Current sample: " + sampleName, TOOL_NAME);
                while ((inputLine = bufferdReader.readLine()) != null) {
                    String line = inputLine.trim();
                    String[] toks = line.split(DELIMITER);
                    if (toks.length == 1) {
                        for (SnpFilter snpFilter : snpFilters) {
                            snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor, minKmers, TOOL_NAME);
                        }
                        samples.add(sampleName);
                        sampleName = line;
                        Reporter.report("[INFO]", "Current sample: " + sampleName, TOOL_NAME);
                    } else if (toks.length != 2) {
                        Reporter.report("[ERROR]", "Unexpected input line: " + line, TOOL_NAME);
                        System.exit(1);
                    }
                    KmerLink kmerLink = map.get(toks[0]);
                    if (kmerLink != null) {
                        boolean setMer = kmerLink.setMer(Short.parseShort(toks[1]));
                        if (!setMer) {
                            System.err.println("mer not set");
                        }
//                    SnpFilter snpFilter = kmerLink.getSnpFilter();
//                    System.err.println(kmerLink.getParentSequence().getId()+"\t"+snpFilter.getSnpPosition()+"\t"+toks[1]);
                    }
                }
                for (SnpFilter snpFilter : snpFilters) {
                    snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor, minKmers, TOOL_NAME);
                }
                samples.add(sampleName);
            } catch (FileNotFoundException ex) {
                Reporter.report("[ERROR]", "File not found: " + kmersFileName, TOOL_NAME);
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
        }
        Reporter.report("[INFO]", "Finished assigning k-mer frequencies to SNPs", TOOL_NAME);
        StringBuilder header = new StringBuilder("Ref:Pos");
        header.append(DELIMITER).append("Phenotype");
        //TODO - risky to just take one!!!
        header.append(DELIMITER).append(parent1);
        header.append(DELIMITER).append(parent2);
        for (String sample : samples) {
            header.append(DELIMITER).append(sample);
        }
        System.out.println(header);

        for (SnpFilter snpFilter : snpFilters) {
            StringBuilder sb = new StringBuilder();
            sb.append(snpFilter.getSequence1().getId()).append("__").append(snpFilter.getSequence2().getId());
            sb.append(":").append(snpFilter.getSnpPosition() + 1).append(DELIMITER);
            sb.append("All").append(DELIMITER);
            sb.append(snpFilter.getBase1()).append(DELIMITER).append(snpFilter.getBase2());
            for (String sample : samples) {
                sb.append(DELIMITER).append(snpFilter.getSnpCall(sample));
                //DEBUG sb.append(" ").append(snpFilter.getSnpCallDetails(sample));
            }
            System.out.println(sb);

//            System.err.println(snpFilter.getSequence1().getId() + DELIMITER + snpFilter.getSequence2().getId()+" CALL: "+snpFilter.callBaseAndResetMers(sampleName, minTotal, minMinor));
//            System.err.println(snpFilter.getMedian1() + " <-- "+Arrays.toString(snpFilter.getMers1()));
//            System.err.println(snpFilter.getMedian2() + " <-- "+Arrays.toString(snpFilter.getMers2()));
        }
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private String intFormat(int value) {
        return NumberFormat.getIntegerInstance().format(value);
    }
}
