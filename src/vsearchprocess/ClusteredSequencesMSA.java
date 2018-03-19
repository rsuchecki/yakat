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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import shared.Reporter;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class ClusteredSequencesMSA {

//    HashMap<String, ArrayList<Sequence>> map;
    private HashMap<String, ClusteredSampleMSA> map;
    private Integer msaAlignmentLength;
    private ArrayList<Snp> intraSnps = new ArrayList<>(0);
    private ArrayList<Snp> interSnps = new ArrayList<>(0);
    private HashMap<String, MsaSeqPair> seqPairs = new HashMap<>();

//    ArrayList<Sequence> list;
    private final String TOOL_NAME;

    public ClusteredSequencesMSA(ArrayList<String> sampleNames, String TOOL_NAME) {
        this.TOOL_NAME = TOOL_NAME;
        map = new HashMap<>(sampleNames.size() * 2);
//        list = new ArrayList<>();
        for (String sampleName : sampleNames) {
            map.put(sampleName, new ClusteredSampleMSA(sampleName));
        }
    }

    public ClusteredSampleMSA getSequences(String sampleName) {
        return map.get(sampleName);
    }

    public ArrayList<MsaSequence> getSequencesList() {
        ArrayList<MsaSequence> sequences = new ArrayList<>();
        for (ClusteredSampleMSA s : map.values()) {
            sequences.addAll(s.getSequences());
        }
        return sequences;
    }

    public boolean hasInterSnps(double minIdentity) {
//        int count = 0;
        for (Snp snp : interSnps) {
            String id1 = snp.getSequence1().getId();
            String id2 = snp.getSequence2().getId();
            String key = id1.compareTo(id2) < 0 ? id1 + id2 : id2 + id1;
            MsaSeqPair pair = seqPairs.get(key);
            if (pair.getMinIdentity() > minIdentity) {
                return true;
//                count++;
            }
        }
//        return count > 0;
        return false;
    }

    public ArrayList<Snp> getIntraSnps() {
        return intraSnps;
    }

//    public ArrayList<Sequence> getSequencesList() {
//        return list;
//    }
    public ArrayList<Snp> getInterSnps() {
        return interSnps;
    }

//
    public int size() {
        int size = 0;
        for (ClusteredSampleMSA s : map.values()) {
            size += s.size();
        }
        return size;
    }

    public Integer getMsaAlignmentLength() {
        return msaAlignmentLength;
    }

    public void addSequence(MsaSequence sequence) {
//        list.add(sequence);
        if (msaAlignmentLength == null) {
            msaAlignmentLength = sequence.getLength();
        } else if (msaAlignmentLength != sequence.getLength()) {
            Reporter.report("[FATAL]", "Sequences in a gievn cluster must be of equal length (with MSA padding where necessary), offending sequence: " + sequence.getId(), TOOL_NAME);
            System.exit(1);
        }

        int added = 0;
        for (String sampleName : map.keySet()) {
            if (sequence.getId().startsWith(sampleName)) {
                if (map.get(sampleName).add(sequence)) {
                    added++;
                }
            }
        }
        if (added > 1) {
            Reporter.report("[FATAL]", "Sequence " + sequence.getId() + " added to more than one subset based on sample/bulk id - possible reason: one id being a substring of another", TOOL_NAME);
            System.exit(1);
        }
        if (added < 1) {
            Reporter.report("[FATAL]", "Sequence " + sequence.getId() + "  not added to a subset based on sample/bulk id - possible reason: id not specified or misspelled", TOOL_NAME);
            System.exit(1);
        }
    }
//
//    /**
//     * Check if all sequences belong to the same sample/bulk
//     *
//     * @return
//     */
//    public boolean isSingleSample() {
//        for (ClusteredSampleMSA sequences : map.values()) {
//            if (sequences.size() == size()) {
//                return true;
//            }
//        }
//        return false;
//    }

    public int getNumClusteredSamples() {
        int samples = 0;
        for (ClusteredSampleMSA c : map.values()) {
            if (c.size() > 0) {
                samples++;
            }
        }
        return samples;
    }

    public boolean mergeSequencesWithinSamples() {
        int merged = 0;
        for (ClusteredSampleMSA sequences : map.values()) {
            if (sequences.mergeSequences()) {
                merged++;
            }
        }
        return merged > 0;
    }

    public void callSNPsBetweenAllSamples(int maxIndelLength, int minIndelDistFromEnds) {
        interSnps = new ArrayList<>(size());
        ClusteredSampleMSA[] samples = map.values().toArray(new ClusteredSampleMSA[map.size()]);
        for (int i = 0; i < samples.length - 1; i++) {
//                System.err.println("["+i+"] "+samples[i].getSampleId()+" VS ");
            for (int j = i + 1; j < samples.length; j++) {
//                System.err.println("\t\t\t["+j+"]"+samples[j].getSampleId());
                callSNPsBetween2Samples(samples[i], samples[j], maxIndelLength, minIndelDistFromEnds);
            }
        }
    }

    public void callSNPsWithinEachSample(int maxIndelLength, int minIndelDistFromEnds) {
        intraSnps = new ArrayList<>(size());
        for (ClusteredSampleMSA s : map.values()) {
            s.callSNPsWithinSample(maxIndelLength, minIndelDistFromEnds);
            ArrayList<Snp> snpsWithin = s.getSnpsWithin();
            if (snpsWithin != null) {
                intraSnps.addAll(snpsWithin);
            }
        }
    }

    private void callSNPsBetween2Samples(ClusteredSampleMSA sample1, ClusteredSampleMSA sample2, int maxIndelLength, int minIndelDistFromEnds) {
        int len = getMsaAlignmentLength();
        for (MsaSequence s1 : sample1.getSequences()) {
            if (s1.getNonTipPaddingLength() > 0) {
                continue;
            }
            String id1 = s1.getId();
            boolean[] padding1 = s1.getPaddingArray(maxIndelLength, minIndelDistFromEnds);
            for (MsaSequence s2 : sample2.getSequences()) {
                String id2 = s2.getId();
                String pairKey = id1.compareTo(id2) < 0 ? id1 + id2 : id2 + id1;
                if (s2.getNonTipPaddingLength() > 0) {
                    continue;
                }
                boolean[] padding2 = s2.getPaddingArray(maxIndelLength, minIndelDistFromEnds);
                for (int i = 0; i < len; i++) {
                    if (s1.getSequenceString().charAt(i) != s2.getSequenceString().charAt(i) && !padding1[i] && !padding2[i]) {
                        Snp snp = new Snp(s1, s2, i);
                        interSnps.add(snp);
                        //COUNT SNPS PER PAIR
                        MsaSeqPair pair = seqPairs.get(pairKey);
                        if (pair == null) {
                            pair = new MsaSeqPair(s1, s2);
                            pair.addSnp(snp);
                            seqPairs.put(pairKey, pair);
                        } else {
                            pair.addSnp(snp);
                        }
                    }
                }
            }
        }
    }

    public ArrayList<Double> getPairwiseIntraIdenities() {
        ArrayList<Double> ids = new ArrayList<>(seqPairs.size());
        for (Map.Entry<String, MsaSeqPair> entry : seqPairs.entrySet()) {
            MsaSeqPair pair = entry.getValue();
            ids.add(pair.getMinIdentity());
        }
        return ids;
    }
//
//    public CharSequence calculateIdentities() {
//        StringBuilder sb = new StringBuilder();
//        for (Map.Entry<String, MsaSeqPair> entry : seqPairs.entrySet()) {
//            MsaSeqPair pair = entry.getValue();
//            sb.append(pair.getS1().getId()).append("\t");
//            sb.append(pair.getS2().getId()).append("\t");            
//            sb.append("MM=").append(pair.getSnps()).append("\t");
//            sb.append("Len1=").append(pair.getS1().getUnpaddedLength()).append("\t");
//            sb.append("Ipd1=").append(pair.getS1().getNonTipPaddingLength()).append("\t");
//            sb.append("Id1=").append(pair.getIdentity1()).append("\t");
//            sb.append("Len2=").append(pair.getS2().getUnpaddedLength()).append("\t");
//            sb.append("Ipd2=").append(pair.getS2().getNonTipPaddingLength()).append("\t");
//            sb.append("Id2=").append(pair.getIdentity2()).append("\t");
//            
//            sb.append(System.lineSeparator());
//        }
//        return sb;
//    }

    public CharSequence getClusterForPrint(int clusterNumber, boolean suppressPadding) {
//        System.out.println("\nCLUSTER: " + clusterLabel);      
        StringBuilder sb = new StringBuilder();
        for (Sequence msaSequence : getSequencesList()) {
            sb.append(">Cluster_").append(clusterNumber).append("_").append(msaSequence.getId());
            sb.append(System.lineSeparator());
            if (suppressPadding) {
                sb.append(msaSequence.getSequenceString().replaceAll("-", "")).append(System.lineSeparator());
            } else {
                sb.append(msaSequence.getSequenceString()).append(System.lineSeparator());
            }
        }
        return sb;
    }

    public void printIntraSnps(int clusterNumber, boolean reverseLex, String DELIMITER,
        String suffix, boolean printSequence) {
        StringBuilder sb = new StringBuilder();
        for (Snp snp : intraSnps) {
            sb.append(snp.getSnpString(clusterNumber, reverseLex, DELIMITER, suffix));
            if (printSequence) {
                sb.append(DELIMITER).append(snp.getSequence1().getSequenceString());
                sb.append(DELIMITER).append(snp.getSequence2().getSequenceString());
            }
            sb.append(System.lineSeparator());
        }
        System.out.print(sb);
    }

    public void printInterSnps(int clusterNumber, boolean reverseLex, String DELIMITER, String suffix,
        double minInterIdentity, boolean printSequence) {
        StringBuilder sb = new StringBuilder();
        for (Snp snp : interSnps) {
            MsaSeqPair pair = seqPairs.get(snp.getSequence1().getId() + snp.getSequence2().getId());
            if (pair.getMinIdentity() > minInterIdentity) {
                sb.append(snp.getSnpString(clusterNumber, reverseLex, DELIMITER, suffix));
                if (printSequence) {
                    sb.append(DELIMITER).append(snp.getSequence1().getSequenceString());
                    sb.append(DELIMITER).append(snp.getSequence2().getSequenceString());
                }
//                sb.append(System.lineSeparator());
//                sb.append(pair.getMinIdentity());
                sb.append(System.lineSeparator());
            }
        }
        System.out.print(sb);
    }

}
