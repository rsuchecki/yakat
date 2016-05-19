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

/**
 *
 * @author rad
 */
public class MsaParserListVariants {

    private final static String DELIMITER = "\t";
//    private char[] getCharsAtPosition(int position) {
//        
//    }
    private int maxIndelLength = 1;
    private int minIndelDistFromEnds = 3;
    private Integer maxSeqsPerCluster;

    public MsaParserListVariants(String fastaFileName) {
        readMSASequencesFromFasta(fastaFileName);
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
