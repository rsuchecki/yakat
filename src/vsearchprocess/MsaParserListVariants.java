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
package vsearchprocess;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import shared.Sequence;
import shared.FastaReader;
import shared.Reporter;
import shared.StdRedirect;

/**
 *
 * @author rad
 */
public class MsaParserListVariants {

    private final static String DELIMITER = "\t";
//    private char[] getCharsAtPosition(int position) {
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int READER_BUFFER_SIZE = 8192;
    private final int WRITER_BUFFER_SIZE = 8192;

//    }
    public MsaParserListVariants(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }
        new StdRedirect(optSet, TOOL_NAME);
        String fileName = (String) optSet.getOpt("clusters-msa").getValueOrDefault();
        readAndProcessMSASequencesFromFasta(fileName, optSet);
//        readMSASequencesFromFasta(optSet.getPositionalOptsList().get(0).getValue());
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Parse MSA output of VSEARCH clustering, call variants within each cluster. "
            + "Variants printed to stdout.");

        //CONSIDER SETTINGS
        //  INPUT LABELS TO DISTINGUISH SAMPLES, PREVENT CALLING SNPS WITHIN SAMPLES
        //  SNP CALLING SETTINGS
        //       -M, --msaFastaSNPsIndels <file>    : Given a FASTA MSA file, list SNPs and indels up to -i, --indelLen bases [REQUIRED for [3]]
        // -x, --maxSeqsPerCluster <int>      : Only output SNPs/indels if up to <int> sequences in the input file, not set by default
        // -i, --maxIndelLen <int>            : Maximum indel length to report, defaults to 1
        // -d, --minIndelDistFromEnds <int>   : Treat indels up to <int> bases from one of the ends as padding and ignore, defaults to 3
        // --sample-id-prefices
        //IF SNPS IN A GIVEN CLUSTER DISPLAY ONE-TO-MANY RELATIONSHIP, CONSIDER MERGING THE "MANY" INTO ONE SEQUENCE
//>*ms5wt_445588
//CATCTTGAGGACGTTGTAGCTATGGTGCGTCACGACCATGAAATTGGCCAGGGTCGGGTACCCGAGGATGACATTGTATGGCAGACGGATGTGAGCGATGTCGAAGTCAATGAGGTCGGTGCGGTAGTTGTCGTGTTGTTCGAAGGTGACGGGGAGGTGGACCTGCCCGATCGGGGTGGTGGAGCCGTCGGTCA
//>ms5mut_154471
//----------------------------------------------GCCAGGGTTGGGTACCCGAGGATGACATTGTATGGCAGACGGATGTGAGCGATGTCGA------------------------------------------------------------------------------------------
//>ms5mut_161305
//----TTGAGGACGTTGTAGCTATGGTGCGTCACGACCATGAAATTGGCCAGGGTTGGGTACC------------------------------------------------------------------------------------------------------------------------------------
//>ms5mut_154274
//-------------------------------------ATGAAATTGGCCAGGGTTGGGTACCCGAGGATGACATTGTATGGCAGACGGATGTG-----------------------------------------------------------------------------------------------------
//>ms5mut_196156
//----------------------------------------------------GTTGGGTACCCGAGGATGACATTGTATGGCAGACGGATGTGAGCGATGTCGAAGTC--------------------------------------------------------------------------------------
//stdin	*ms5wt_445588	ms5mut_154471	55	C	T
//stdin	*ms5wt_445588	ms5mut_161305	55	C	T
//stdin	*ms5wt_445588	ms5mut_154274	55	C	T
//stdin	*ms5wt_445588	ms5mut_196156	55	C	T
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt(null, "sample-ids", "Space separated sample identifiers which form the prefices of the input FASTA identifiers")
            .setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE).setRequired(true));
        optSet.addOpt(new Opt(null, "clusters-msa", "The vsearch cluster msaout file, alternatively use stdin", 1));
//        optSet.addOpt(new Opt(null, "original-fasta", "The original FASTA file gven to vsearch, if specified it will be used to extract the input sequences' lengths", 1));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Base calling settings]");
//        optSet.addOpt(new Opt('a', "min-coverage-per-allele", "Minimum coverage required for an allele to be considered in a locus call", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt('A', "max-percent-error-allele", "Percentage of coverage up to which an allele is regarded to be an error", 1).setMinValue(0.0).setDefaultValue(1.0));
//        optSet.addOpt(new Opt('L', "max-percent-error-locus", "Percentage coverage of alternative alleles up to which a locus is reported", 1).setMinValue(0.0).setDefaultValue(1.0));
//        optSet.addOpt(new Opt(null, "min-minor-major-ratio", "[TODO] Minimum fraction of a minor allele bases required to call a heterozygous base", 1).setMinValue(0.0).setMaxValue(0.5));
//        optSet.addOpt(new Opt(null, "zero-reads-char", "A character denoting zero reads at a postion for a given sample", 1).setDefaultValue('.'));
//        optSet.addOpt(new Opt(null, "ambiguous-call-char", "A character indicating uncertain call (e.g. due to low coverage or unclear zygosity at locus)", 1).setDefaultValue('?'));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Input clusters processing settings]");
        optSet.addOpt(new Opt(null, "min-samples-clustered", "Minimum number of samples in an input  cluster", 1).setMinValue(1).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "min-seqs-clustered-in", "Minimum number of sequences in an input cluster", 1).setMinValue(2).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "max-seqs-clustered-in", "Maximum number of sequences in an input cluster", 1).setMinValue(2).setDefaultValue(1000));
        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Variant calling and reporting]");
//        optSet.addOpt(new Opt(null, "min-sequences-per-cluster", "Minimum number of sequences required for a cluster to be considered ", 1).setMinValue(2).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "max-seqs-clustered-out", "Maximum number of sequences in an output cluster", 1).setMinValue(2).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "min-inter-identity", "Minimum inter sample identity ", 1).setMinValue(0.0).setDefaultValue(0.95));
        optSet.addOpt(new Opt(null, "max-inter-snps", "SNPs will be reported if at most <arg> inter-sample SNPs are called in a cluster", 1).setMinValue(0).setDefaultValue(2));
        optSet.addOpt(new Opt(null, "max-intra-snps", "SNPs will be reported if at most <arg> intra-sample SNPs are called in a cluster", 1).setMinValue(0).setDefaultValue(1));
        optSet.addOpt(new Opt(null, "max-indel-length", "Maximum length of an indel ", 1).setMinValue(0).setDefaultValue(1));
        optSet.addOpt(new Opt(null, "min-indel-distance", "Treat '-' positions less than <arg> bases from sequence ends as padding not indels", 1).setMinValue(1).setDefaultValue(3));
        optSet.addOpt(new Opt(null, "supress-inter-snps", "Do not report inter-sample SNPs"));
        optSet.addOpt(new Opt(null, "supress-intra-snps", "Do not report intra-sample SNPs"));
        optSet.addOpt(new Opt(null, "reverse-lex-order", "Variants will be reported in reverse lexicographical order of the sample/bulk identifiers"));

        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Runtime and output settings]");
//        String threadsOrderNote = "Note that in multi-threaded mode the output lines order need not reflect the input order";
//        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(1, threadsOrderNote));
//        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 1024, 128, 32768));
//        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",64, 1, 256));
////        optSet.addOpt(new Opt('u', "out-buffer-size", "Size of buffers put on out-queue ", 1024, 128, 32768));
////        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writing-out",64, 1, 256));
        optSet.addOpt(new Opt(null, "out-clusters-msa", "Output clustered sequences (for which SNPs were called) to <arg> MSA/FASTA file", 1));
        optSet.addOpt(new Opt(null, "out-unclustered-fasta", "Output unclustered sequences to <arg> FASTA file", 1));
        optSet.addOpt(new Opt(null, "out-unclustered-min-len", "Minimum length required to output an unclustered sequence", 1).setMinValue(1).setDefaultValue(100));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
////        String headerNote = "Can be useful for external parallization (print header once)";
////        optSet.addOpt(new Opt('H', "header-only", "Print header and exit").addFootnote(1, TOOL_NAME));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[A little bit of help]");
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
//        optSet.addOpt(new Opt('I', "iupac-codes-table", "Print the table of IUPAC nucleotide codes and exit"));
//        optSet.addOpt(new Opt('D', "additional-codes-table", "Print the table of additional codes/symbols used by this program"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    public void readAndProcessMSASequencesFromFasta(String fileName, OptSet optSet) {
        ArrayList<String> SAMPLE_NAMES = (ArrayList<String>) optSet.getOpt("sample-ids").getValues();

//        String originalFastaFileName = (String)optSet.getOpt("original-fasta").getValueOrDefault();
//        HashMap<String, Integer> lengthsMap;
//        if(originalFastaFileName != null) {
//            lengthsMap = FastaReader.hashMapOfSequenceLengthsFromFasta(originalFastaFileName, null); 
//        } else {
//            lengthsMap = new HashMap<>(0);
//        }
//        ArrayList<Sequence> sequencesList = new ArrayList<>();
        StringBuilder sb = new StringBuilder("ClusterId");
        sb.append(DELIMITER);
        sb.append("AlnLen");
        sb.append(DELIMITER);
        sb.append("Id1");
        sb.append(DELIMITER);
        sb.append("Len1");
        sb.append(DELIMITER);
        sb.append("Id2");
        sb.append(DELIMITER);
        sb.append("Len2");
        sb.append(DELIMITER);
        sb.append("Pos");
        sb.append(DELIMITER);
        sb.append("Base1");
        sb.append(DELIMITER);
        sb.append("Base2");
        sb.append(DELIMITER);
        sb.append("Comments");
        sb.append(DELIMITER);
        sb.append("Seq1");
        sb.append(DELIMITER);
        sb.append("Seq2");
        System.out.println(sb);

        ClusteredSequencesMSA clusteredSeqs = new ClusteredSequencesMSA(SAMPLE_NAMES, TOOL_NAME);
        BufferedReader bufferdReader = null;

        BufferedWriter unclusteredFastaOut = null;
        BufferedWriter clustersFastaOut = null;
        String unclusteredOutFile = (String) optSet.getOpt("out-unclustered-fasta").getValueOrDefault();
        String clustersOutFile = (String) optSet.getOpt("out-clusters-msa").getValueOrDefault();
        int unclusteredOutMinLength = (int) optSet.getOpt("out-unclustered-min-len").getValueOrDefault();
        try {
            if (unclusteredOutFile != null) {
                unclusteredFastaOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(unclusteredOutFile))), WRITER_BUFFER_SIZE);
            }
            if (clustersOutFile != null) {
                clustersFastaOut = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(clustersOutFile))), WRITER_BUFFER_SIZE);
            }
            String inputLine;
            if (fileName == null) {
//                Reporter.report("[INFO]", "Input file(s) not specified, reading from stdin ", TOOL_NAME);
                bufferdReader = new BufferedReader(new InputStreamReader(System.in), READER_BUFFER_SIZE);
            } else if (fileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(fileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(fileName)), READER_BUFFER_SIZE);
            }
            String id = "";
            StringBuilder seqBuilder = new StringBuilder();
            int clusterNumber = 0;
            while ((inputLine = bufferdReader.readLine()) != null) {
                String line = inputLine.trim();
                if (line.startsWith(">")) {
                    if (seqBuilder.length() > 0) {   //SEQUENCE STORED EARLIER                        
//                        sequencesList.add(new Sequence(id, seqBuilder.toString()));
                        if (!id.startsWith("consensus")) {
                            clusteredSeqs.addSequence(new MsaSequence(id, seqBuilder.toString()));
//                            clusteredSeqs.addSequence(new MsaSequence(id, seqBuilder.toString(), lengthsMap.get(id)));
                        }
                        seqBuilder = new StringBuilder();
                    }
                    id = line.substring(1); //store current id, get rid of ">"
                    if (line.startsWith(">*")) { //INDICATING NEW CLUSTER                        
                        id = id.substring(1); //get rid of "*"
                        //INIT NEW 
//                        sequencesList = new ArrayList<>();
                        clusteredSeqs = new ClusteredSequencesMSA(SAMPLE_NAMES, TOOL_NAME);
                    } else if (line.equals(">consensus")) {
                        //PROCESS PREVIOUS CLUSER 
//                        if (sequencesList.size() > 1) {
                        clusterNumber = processCluster(clusteredSeqs, optSet, clusterNumber, clustersFastaOut);

                        if (unclusteredFastaOut != null && clusteredSeqs.size() == 1) {
                            MsaSequence seq = clusteredSeqs.getSequencesList().get(0);
                            if (seq.getUnpaddedLength() >= unclusteredOutMinLength) {
                                unclusteredFastaOut.write(clusteredSeqs.getSequencesList().get(0).getFasta(true).toString());
                                unclusteredFastaOut.newLine();
                            }
                        }
//                        }
                        //SKIP the consensus
//                        continue;
                    }

                } else {
                    seqBuilder.append(line);
                }
            }
            if (unclusteredFastaOut != null) {
                unclusteredFastaOut.flush();
                unclusteredFastaOut.close();
            }
            if (clustersFastaOut != null) {
                clustersFastaOut.flush();
                clustersFastaOut.close();
            }

            if (id.isEmpty()) {
                System.err.println("Error reading FASTA [" + fileName + "]. No identifier lines? Terminating... ");
                System.exit(1);
            } else {
//                System.err.println("Can ignore this if don't care about the consensus sequence");
//                String identifierString[] = id.split(" ");
//                String key = identifierString[0];
//                    sequencesList.add(new Sequence(id, seqBuilder.toString()));
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + fileName, TOOL_NAME);
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

    private int processCluster(ClusteredSequencesMSA clusteredSeqs, OptSet optSet, int clusterNumber, BufferedWriter clustersFastaOut) throws IOException {

        int minSamplesClustered = (int) optSet.getOpt("min-samples-clustered").getValueOrDefault();
        int minSeqsClusteredIn = (int) optSet.getOpt("min-seqs-clustered-in").getValueOrDefault();
        int maxSeqsClusteredIn = (int) optSet.getOpt("max-seqs-clustered-in").getValueOrDefault();
        
        int maxSeqsClusteredOut = (int) optSet.getOpt("max-seqs-clustered-out").getValueOrDefault();
        int maxIntraSnps = (int) optSet.getOpt("max-intra-snps").getValueOrDefault();
        int maxInterSnps = (int) optSet.getOpt("max-inter-snps").getValueOrDefault();

        int maxIndelLength = (int) optSet.getOpt("max-indel-length").getValueOrDefault();
        int minIndelDistFromEnds = (int) optSet.getOpt("min-indel-distance").getValueOrDefault();
        
        double minInterIdentity = (double) optSet.getOpt("min-inter-identity").getValueOrDefault();

        boolean reverseLex = optSet.getOpt("reverse-lex-order").getOptFlag();
        boolean supressIntra = optSet.getOpt("supress-intra-snps").getOptFlag();
        boolean supressInter = optSet.getOpt("supress-inter-snps").getOptFlag();

        boolean appendSequencesToSnpList = true;
        
        if (clusteredSeqs.size() >= minSeqsClusteredIn && clusteredSeqs.size() <= maxSeqsClusteredIn && clusteredSeqs.getNumClusteredSamples() >= minSamplesClustered) {
            //CALL WITHIN EACH SAMPLE
            clusteredSeqs.callSNPsWithinEachSample(maxIndelLength, minIndelDistFromEnds);
            String suffix = null;
            boolean hasIntra = false;
            if (!clusteredSeqs.getIntraSnps().isEmpty()) {
                suffix = "HAS_INTRA";
                hasIntra = true;
            }
//            String clusterString = "ALL";
            //MERGE NON-CONFLICTING SEQUENCES WITHIN EACH SAMPLE
            if (clusteredSeqs.mergeSequencesWithinSamples()) {
//                clusterString = "MERGED";
            }
            
            //CALL BETWEEN SAMPLES
            clusteredSeqs.callSNPsBetweenAllSamples(maxIndelLength, minIndelDistFromEnds);
//            boolean hasInter = false;
//            if (!clusteredSeqs.getInterSnps().isEmpty()) {
//                hasInter = true;
//            }
            boolean hasInter = clusteredSeqs.hasInterSnps(minInterIdentity);
//            ArrayList<Double> pairwiseIntraIdenities = clusteredSeqs.getPairwiseIntraIdenities();
//            try {
//            Double minIdentity = Collections.min(pairwiseIntraIdenities);
//            System.out.println(minIdentity);
//            } catch (NoSuchElementException e) {
//                int x =0;
//            }
            
            int intra = clusteredSeqs.getIntraSnps().size();
            int inter = clusteredSeqs.getInterSnps().size();
            //PRINT            
            if (intra <= maxIntraSnps &&  inter <= maxInterSnps && clusteredSeqs.size() <= maxSeqsClusteredOut) {
//                clusteredSeqs.printCluster((clusterNumber) + " " + clusterString, maxIndelLength);
                if (clustersFastaOut != null && ((hasIntra && !supressIntra) || (hasInter && !supressInter))) {
                    clustersFastaOut.write(clusteredSeqs.getClusterForPrint(++clusterNumber).toString());
//                    clustersFastaOut.newLine();
                }
                if (!supressIntra) {
                    clusteredSeqs.printIntraSnps(clusterNumber, reverseLex, DELIMITER, "INTRA", appendSequencesToSnpList);
                }
                if (!supressInter) {
                    clusteredSeqs.printInterSnps(clusterNumber, reverseLex, DELIMITER, suffix, minInterIdentity, appendSequencesToSnpList);
                }
            }
        }
        return clusterNumber;
    }

//    private class SequenceLengthComparator implements Comparator<Sequence> {
//
//        @Override
//        public int compare(Sequence sequence, Sequence anotherSequence) {
//            return sequence.getLength() - anotherSequence.getLength();
//        }
//    }
}
