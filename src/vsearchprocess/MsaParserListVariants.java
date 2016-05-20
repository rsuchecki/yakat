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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Set;
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

//    }
    private int maxIndelLength = 1;
    private int minIndelDistFromEnds = 3;
    private Integer maxSeqsPerCluster;

    public MsaParserListVariants(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        new StdRedirect(optSet, TOOL_NAME);
        readMSASequencesFromFasta(optSet.getPositionalOptsList().get(0).getValue());
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Parse MSA output of VSEARCH clustering, call variants within each cluster. ");
        
        //CONSIDER SETTINGS
        //  INPUT LABELS TO DISTINGUISH SAMPLES, PREVENT CALLING SNPS WITHIN SAMPLES
        //  SNP CALLING SETTINGS
                //       -M, --msaFastaSNPsIndels <file>    : Given a FASTA MSA file, list SNPs and indels up to -i, --indelLen bases [REQUIRED for [3]]
                // -x, --maxSeqsPerCluster <int>      : Only output SNPs/indels if up to <int> sequences in the input file, not set by default
                // -i, --maxIndelLen <int>            : Maximum indel length to report, defaults to 1
                // -d, --minIndelDistFromEnds <int>   : Treat indels up to <int> bases from one of the ends as padding and ignore, defaults to 3

        
        
        
        //INPUT
//        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('n', "sample-names", "space separated sample names, order must correspond to mpileup input").setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE));
//        optSet.addOpt(new Opt('p', "pileup-file", "Input (m)pileup file, alternatively use stdin", 1));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Base calling settings]");
//        optSet.addOpt(new Opt('a', "min-coverage-per-allele", "Minimum coverage required for an allele to be considered in a locus call", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt('c', "min-coverage-per-locus", "Minimum coverage required for a locus to be considered ", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt('C', "max-coverage-per-locus", "Maximum coverage allowed for a locus to be considered", 1).setMinValue(1).setDefaultValue(1000));
//        optSet.addOpt(new Opt('A', "max-percent-error-allele", "Percentage of coverage up to which an allele is regarded to be an error", 1).setMinValue(0.0).setDefaultValue(1.0));
//        optSet.addOpt(new Opt('L', "max-percent-error-locus", "Percentage coverage of alternative alleles up to which a locus is reported", 1).setMinValue(0.0).setDefaultValue(1.0));
//        optSet.addOpt(new Opt(null, "min-minor-major-ratio", "[TODO] Minimum fraction of a minor allele bases required to call a heterozygous base", 1).setMinValue(0.0).setMaxValue(0.5));
//        optSet.addOpt(new Opt(null, "zero-reads-char", "A character denoting zero reads at a postion for a given sample", 1).setDefaultValue('.'));
//        optSet.addOpt(new Opt(null, "ambiguous-call-char", "A character indicating uncertain call (e.g. due to low coverage or unclear zygosity at locus)", 1).setDefaultValue('?'));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[Record reporting settings]");
//        optSet.addOpt(new Opt('s', "min-samples-within-coverage", "Minimum samples within coverage thresholds required to produce ouput", 1).setMinValue(0).setDefaultValue(2));
//        
//        optSet.addOpt(new Opt('W', "all-within-thresholds", "Print all loci within given thresholds even if no alternative alleles called"));
//        optSet.addOpt(new Opt(null, "min-samples-called", "[TODO] Minimum samples for which the base was called", 1).setMinValue(1).setDefaultValue(2));
//        optSet.addOpt(new Opt(null, "max-samples-called", "[TODO] Maximum samples for which the base was called", 1).setMinValue(1));
//        optSet.addOpt(new Opt(null, "min-samples-zero-coverage", "[TODO] May be useful for presence-absence analyses", 1).setMinValue(0).setDefaultValue(0));
//        optSet.addOpt(new Opt(null, "max-samples-zero-coverage", "[TODO] May be useful for presence-absence analyses", 1).setMinValue(0));
//        
//        optSet.addOpt(new Opt(null, "min-snps-to-reference", "[TODO] Report a position if at least <arg> samples have a SNP to reference ",1 ).setMinValue(0).setDefaultValue(0));
//        optSet.addOpt(new Opt(null, "min-calls-uncertain", "Minimum samples for which the base was not called due to uncertainty", 1).setMinValue(0).setDefaultValue(0));
//        optSet.addOpt(new Opt(null, "max-calls-uncertain", "Maximum samples for which the base was not called due to uncertainty", 1).setMinValue(0).setDefaultValue(65535));
//        optSet.addOpt(new Opt(null, "min-calls-het", "Report a position if at least <arg> calls are heterozygous",1 ).setMinValue(0).setDefaultValue(0));
//        optSet.addOpt(new Opt(null, "max-calls-het", "Report a position if at most <arg> calls are heterozygous",1 ).setMinValue(0).setDefaultValue(65535));
////        optSet.addOpt(new Opt('z', "min-missing-samples", "Minimum samples with zero coverage", 1).setMinValue(0).setDefaultValue(0));
////        optSet.addOpt(new Opt('u', "max-uncalled-samples", "Maximum samples for which the base was not called", 1).setMinValue(0).setDefaultValue(0));
//        optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel("[Runtime settings]");
//        String threadsOrderNote = "Note that in multi-threaded mode the output lines order need not reflect the input order";
//        optSet.addOpt(new Opt('t', "threads", "Max number of threads to be used", 1).setMinValue(1).setDefaultValue(1).setMaxValue(Runtime.getRuntime().availableProcessors()).addFootnote(1, threadsOrderNote));
//        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 1024, 128, 32768));
//        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",64, 1, 256));
////        optSet.addOpt(new Opt('u', "out-buffer-size", "Size of buffers put on out-queue ", 1024, 128, 32768));
////        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writing-out",64, 1, 256));
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
////        String headerNote = "Can be useful for external parallization (print header once)";
////        optSet.addOpt(new Opt('H', "header-only", "Print header and exit").addFootnote(1, TOOL_NAME));
//        optSet.incrementLisitngGroup();
//        optSet.setListingGroupLabel("[A little bit of help]");
//        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
//        optSet.addOpt(new Opt('I', "iupac-codes-table", "Print the table of IUPAC nucleotide codes and exit"));
//        optSet.addOpt(new Opt('D', "additional-codes-table", "Print the table of additional codes/symbols used by this program"));
        boolean positionalArgumentRequired = true;
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }
    
    public void readMSASequencesFromFasta(String fastaFileName) {
        ArrayList<Sequence> sequencesList = new ArrayList<>();
        File newFile = new File(fastaFileName);
        BufferedReader bufferdReader = null;
        try {
            String inputLine;
            bufferdReader = new BufferedReader(new FileReader(newFile));
            String id = "";
            StringBuilder seqBuilder = new StringBuilder();
            while ((inputLine = bufferdReader.readLine()) != null) {
                String line = inputLine.trim();
                if (line.startsWith(">")) {
                    if (seqBuilder.length() > 0) {   //SEQUENCE STORED EARLIER                        
                        sequencesList.add(new Sequence(id, seqBuilder.toString()));
                        seqBuilder = new StringBuilder();
                    }
                    id = line.substring(1); //store current id, get rid of ">"
                    if (line.startsWith(">*")) { //INDICATING NE CLUSTER                        
                        //INIT NEW 
                        sequencesList = new ArrayList<>();
                    } else if(line.equals(">consensus")) {
                        //PROCESS PREVIOUS CLUSER 
                        if(sequencesList.size() > 1) {
                            parseMSA(sequencesList);
                        }
                        //SKIP the consensus
//                        continue;
                    }

                } else {
                    seqBuilder.append(line);
                }
            }
            if (id.isEmpty()) {
                System.err.println("Error reading FASTA [" + fastaFileName + "]. No identifier lines? Terminating... ");
                System.exit(1);
            } else {
                System.err.println("Can ignore this if dont' care about the consensus sequence");
//                String identifierString[] = id.split(" ");
//                String key = identifierString[0];
//                    sequencesList.add(new Sequence(id, seqBuilder.toString()));
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + newFile.getName(), "PUT-NAME-HERE");
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

//    public void parseMSA(ArrayList<Sequence> msaSequences, int maxIndelLength, int minIndelDistFromEnds, Integer maxSeqsPerCluster) {
    public void parseMSA(ArrayList<Sequence> msaSequences) {
        System.out.println("CLUSTER");
        for (Sequence msaSequence : msaSequences) {
            System.out.println(">"+msaSequence.getId());
            System.out.println(msaSequence.getSequenceString());
        }
        
        
//        HashMap<String, ArrayList<Integer>> idToIndelPostion = new HashMap<>();
//        ArrayList<Sequence> msaSequences = FastaReader.listOfSequencesFromFasta(filenameMsaOfACluster, null);        
//        if(msaSequences.isEmpty()) {
//            System.err.println("Error. No FASTA records read from expected MSA file: "+filenameMsaOfACluster);
//            System.exit(1);
//        }
        int numSeqs = msaSequences.size();
        if (maxSeqsPerCluster == null || numSeqs <= maxSeqsPerCluster) { //if set, ignore clusters with more than maxSeqsPerCluster elements 
            Collections.sort(msaSequences, new SequenceLengthComparator());

            int alnLen = msaSequences.get(0).getLength();

            //record positions regarded as padding (rather than indels)
            ArrayList<boolean[]> paddingsList = new ArrayList<>(numSeqs);
            for (Sequence sequence : msaSequences) {
                boolean[] padding = new boolean[alnLen];
//            System.err.println(sequence.getSequenceString());
                for (int j = 0; j < alnLen; j++) {
                    if (sequence.getSequenceString().charAt(j) == '-') {
                        padding[j] = true;
//                    System.err.println("Padding pos = " + (j + 1) + " in seq " + sequence.getIdentifierString());
                    }
                }
                int padLength = 0;
                for (int positionInSeq = 0; positionInSeq < padding.length + 1; positionInSeq++) {
                    boolean processLast = false;
                    if (positionInSeq == padding.length) {  //Getlast position
                        processLast = true;
                    }

                    if (!processLast && padding[positionInSeq]) {
                        padLength++;
                    } else {

                        //IF padding strech shorter than maxIndelLength THEN it may be an indel and not padding
                        if (padLength > 0 && padLength <= maxIndelLength) {
                            boolean headPadding = false;
                            for (int positionInPadStretch = padLength; positionInPadStretch > 0; positionInPadStretch--) {
                                //IGNORE indels adjacent to SEQ ENDS, that is, treat them as alignment padding
                                if (positionInSeq - positionInPadStretch < minIndelDistFromEnds || headPadding || positionInSeq > alnLen - minIndelDistFromEnds) {
                                    headPadding = true;
//                                System.err.println("Not unpadding pos = " + (positionInSeq - positionInPadStretch + 1) + " in seq " + sequence.getIdentifierString());
                                    break; //IF WE CLASS THIS POSITION AS PADDING, THE REST OF THE STRETCH IS PADDING AS WELL                                
                                } else {
                                    padding[positionInSeq - positionInPadStretch] = false;
                                    headPadding = false;
//                                System.err.println("Unpadding pos = " + (positionInSeq - positionInPadStretch + 1) + " in seq " + sequence.getIdentifierString());
                                }
                            }
                        }
                        padLength = 0;
                    }
                }
                paddingsList.add(padding);
            }

            //FOR EACH POSITION
            for (int position = 0; position < alnLen; position++) {
                char[] bases = new char[numSeqs];
                //RECORD BASES FROM ALL SEQS
                for (int i = 0; i < numSeqs; i++) {
                    bases[i] = msaSequences.get(i).getSequenceString().charAt(position);
                }
                //FOR RECORDED BASE OF A SEQUENCE, COMPARE WITH BASES IN OTHER SEQS
                for (int seq = 0; seq < bases.length - 1; seq++) {
                    char base = bases[seq];
                    boolean padding = paddingsList.get(seq)[position];
                    for (int anotherSeq = seq + 1; anotherSeq < bases.length; anotherSeq++) {
                        boolean anotherSeqPadding = paddingsList.get(anotherSeq)[position];
                        if (base != bases[anotherSeq] && !padding && !anotherSeqPadding) {
                            String source;
//                            if (filenameMsaOfACluster == null) {
                            source = "stdin";
//                            } else {
//                                source = filenameMsaOfACluster.replaceFirst(".*/", "");
//                            }
                            StringBuilder sb = new StringBuilder(source);

                            String id1 = msaSequences.get(seq).getId();
                            String id2 = msaSequences.get(anotherSeq).getId();

                            if (id1.compareTo(id2) < 0) {
                                sb.append(DELIMITER).append(id1);
                                sb.append(DELIMITER).append(id2);
                                sb.append(DELIMITER).append(position + 1).append(DELIMITER).append(base).append(DELIMITER).append(bases[anotherSeq]);
                            } else {
                                sb.append(DELIMITER).append(id2);
                                sb.append(DELIMITER).append(id1);
                                sb.append(DELIMITER).append(position + 1).append(DELIMITER).append(bases[anotherSeq]).append(DELIMITER).append(base);

                            }

                            System.out.println(sb.toString());
                        }

                    }

                }
            }
        }

    }

//    /**
//     * Corrects indel positions for offset introduced by MSA and where
//     * applicable reverse complementing before MSA so that the corrected
//     * position matches the original non rc'd sequence Suffix '_rc' is removed
//     * if present in the seq identifier
//     *
//     * @param indelPosition
//     * @param sequence
//     * @return
//     */
//    private int getCorrectedIndelPosition(Integer indelPosition, Sequence sequence) {
//        String seqUpToIndel;
//        //WHAT IF SEQ RC'd WHEN PARSING CD-HIT OUTPUT?
//        //1. 
//        //2. REMOVE rc suffix from identifier
//        //3. RC INDEl POSITION
//        //4. PROCEED WITH COORDNIATE CORRECTION 
//        if (sequence.getSequenceString().charAt(indelPosition - 1) != '-') {
//            indelPosition--; //offset fix where a '-' was not introduced in the sequence
//        }
//        if (sequence.getIdentifierString().endsWith("_rc")) {
//            seqUpToIndel = sequence.getSequenceString().substring(indelPosition);
////            indelPosition = seqUpToIndel.length();
//            sequence.setIdentifierString(sequence.getIdentifierString().replaceAll("_rc", "")); //SMALL _rc  ONLY, THE CAPITAL _RC HAS NOTHING TO DO WITH OUR CH-HIT OR MSA!
//        } else {
//            seqUpToIndel = sequence.getSubsequence(0, indelPosition);
//        }
//        String ungappedUpToIndel = seqUpToIndel.replaceAll("-", "");
//        return ungappedUpToIndel.length();
//    }
    private class SequenceLengthComparator implements Comparator<Sequence> {

        @Override
        public int compare(Sequence sequence, Sequence anotherSequence) {
            return sequence.getLength() - anotherSequence.getLength();
        }
    }
}
