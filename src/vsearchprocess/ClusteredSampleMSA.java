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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

/**
 *
 * @author rad
 */
public class ClusteredSampleMSA {

    private final ArrayList<MsaSequence> sequences;
    private final String sampleId;
    private ArrayList<Snp> snpsWithin;

    public ClusteredSampleMSA(String sample) {
        sequences = new ArrayList<>();
        this.sampleId = sample;
    }

    public int size() {
        return sequences.size();
    }

    public String getSampleId() {
        return sampleId;
    }

    public boolean add(MsaSequence s) {
        return sequences.add(s);
    }

    public ArrayList<MsaSequence> getSequences() {
        return sequences;
    }

    public ArrayList<Snp> getSnpsWithin() {
        return snpsWithin;
    }

    public int getSnpsWithinCount() {
        return snpsWithin == null ? 0 : snpsWithin.size();
    }

    public boolean mergeSequences() {
        boolean merged = false;
        int numSeqs = getSequences().size();
        if (numSeqs > 1) {
//                ArrayList<boolean[]> paddingList = getPaddingList(sequences, maxIndelLength, minIndelDistFromEnds);
            boolean mismatch = false;
            int alnLen = getSequences().get(0).getLength();
            StringBuilder consensus = new StringBuilder();
            //FOR EACH POSITION
//            compare:
//            for (int position = 0; position < alnLen; position++) {
//                char[] bases = new char[numSeqs];
//
//                //RECORD BASES FROM ALL SEQS
//                for (int i = 0; i < numSeqs; i++) {
//                    bases[i] = getSequences().get(i).getSequenceString().charAt(position);
//                }
//                consensus.append(getBase(bases));
//                //FOR RECORDED BASE OF A SEQUENCE, COMPARE WITH BASES IN OTHER SEQS                    
//                for (int seq = 0; seq < bases.length - 1; seq++) {
//                    char base = bases[seq];
////                        boolean padding = paddingList.get(seq)[position];
//                    for (int anotherSeq = seq + 1; anotherSeq < bases.length; anotherSeq++) {
////                            boolean anotherSeqPadding = paddingList.get(anotherSeq)[position];
//                        if (base != bases[anotherSeq] && base != '-' && bases[anotherSeq] != '-') {
////                            if (base != bases[anotherSeq] && !padding && !anotherSeqPadding) {
//                            mismatch = true;
//                            break compare;
//                        }
//                    }
//                }
//            }
//            consensus = haveMismatches(false);
            //WE MIGHT HAVE MISSED MISMATCHES DUE TO INTERNAL PADDING 
//            consensus = haveMismatches(true) == null ? null : consensus;
            //CREATE NEW MsaSequence OBJECTS AS PADDING ARRAY IS PERSISTANT AND FOR THIS PURPOSE WE ONLY TREAT TERMINAL DASHES AS PADDING
            ArrayList<MsaSequence> seqs = new ArrayList<>(numSeqs);
            for (MsaSequence sequence : getSequences()) {
                MsaSequence msaSequence = new MsaSequence(sequence.getId(), sequence.getSequenceString());
                seqs.add(msaSequence);
//                boolean[] paddingArray = msaSequence.getPaddingArray(Integer.MAX_VALUE, 1);
//                System.err.println(Arrays.toString(paddingArray));
//                int x = 0;
            }

            compare:
            for (int position = 0; position < alnLen; position++) {
                char bases[] = new char[numSeqs];
                boolean padding[] = new boolean[numSeqs];

                //RECORD BASES FROM ALL SEQS
                for (int i = 0; i < numSeqs; i++) {
                    bases[i] = seqs.get(i).getSequenceString().charAt(position);                    
                    padding[i] = seqs.get(1).getPaddingArray(Integer.MAX_VALUE, 0)[position];
                }
                consensus.append(getBase(bases));
                //FOR RECORDED BASE OF A SEQUENCE, COMPARE WITH BASES IN OTHER SEQS                    
                for (int seq = 0; seq < bases.length - 1; seq++) {
                    char base = bases[seq];
                    boolean pad = padding[seq];
//                        boolean padding = paddingList.get(seq)[position];
                    for (int anotherSeq = seq + 1; anotherSeq < bases.length; anotherSeq++) {
//                            boolean anotherSeqPadding = paddingList.get(anotherSeq)[position];
//                    if (base != bases[anotherSeq] && base != '-' && bases[anotherSeq] != '-') {
                        //MISMATCH & NEIGHTER OF POSITION IS TERMINAL PADDING
                        if (base != bases[anotherSeq] && !pad && !padding[anotherSeq]) {
//                            if (base != bases[anotherSeq] && !padding && !anotherSeqPadding) {
                            mismatch = true;
//                            break compare;
                        }
                    }
                }
            }

            if (!mismatch) {            
                StringBuilder id = new StringBuilder(getSampleId());
                for (MsaSequence sequence : getSequences()) {
                    id.append(sequence.getId().replaceFirst(getSampleId(), ""));
//                    list.remove(sequence); //REMOVE FROM THE GLOBAL LIST
                }

                MsaSequence msaSequence = new MsaSequence(id.toString(), consensus.toString());

                //IF PADDING WITHIN MERGED SEQUENCE -> DO NOT OVERWRITE THE ORIGINAL SEQS
                if (msaSequence.getNonTipPaddingLength() > 0) {
//                    System.out.println("PADDING: "+Arrays.toString(split));                    
                } else {

                    getSequences().clear();
                    getSequences().add(msaSequence);
//                ArrayList<Sequence> consSequenceList = new ArrayList<>(1);
//                consSequenceList.add(sequence);
//                map.put(sampleId, consSequenceList); //REPLACE IN MAP
//                list.add(sequence); //ADD TO GLOBAL LIST
                    merged = true;
                }
            }
        }
        return merged;
    }

    private char getBase(char[] bases) {
        for (char base : bases) {
            if (base != '-') {
                return base;
            }
        }
        return '-';
    }

    private StringBuilder haveMismatches() {
        int numSeqs = getSequences().size();
        int alnLen = getSequences().get(0).getLength();
        StringBuilder consensus = new StringBuilder();
        ArrayList<MsaSequence> seqs = new ArrayList<>(numSeqs);
        for (MsaSequence sequence : getSequences()) {
            MsaSequence msaSequence = new MsaSequence(sequence.getId(), sequence.getSequenceString());
            boolean[] paddingArray = msaSequence.getPaddingArray(Integer.MAX_VALUE, 1);
            System.err.println(Arrays.toString(paddingArray));
            int x = 0;
        }

        for (int position = 0; position < alnLen; position++) {
            char bases[] = new char[numSeqs];
            boolean padding[] = new boolean[numSeqs];

            //RECORD BASES FROM ALL SEQS
            for (int i = 0; i < numSeqs; i++) {
                bases[i] = seqs.get(i).getSequenceString().charAt(position);
                padding[i] = seqs.get(1).getPaddingArray(Integer.MAX_VALUE, 1)[position];
            }
            consensus.append(getBase(bases));
            //FOR RECORDED BASE OF A SEQUENCE, COMPARE WITH BASES IN OTHER SEQS                    
            for (int seq = 0; seq < bases.length - 1; seq++) {
                char base = bases[seq];
                boolean pad = padding[seq];
//                        boolean padding = paddingList.get(seq)[position];
                for (int anotherSeq = seq + 1; anotherSeq < bases.length; anotherSeq++) {
//                            boolean anotherSeqPadding = paddingList.get(anotherSeq)[position];
//                    if (base != bases[anotherSeq] && base != '-' && bases[anotherSeq] != '-') {
                    if (base != bases[anotherSeq] && !pad && !padding[anotherSeq]) {
//                            if (base != bases[anotherSeq] && !padding && !anotherSeqPadding) {
//                            mismatch = true;
                        return null;
//                            break compare;
                    }
                }
            }
        }
        return consensus;
    }

    public void callSNPsWithinSample(int maxIndelLength, int minIndelDistFromEnds) {
        if (size() > 1) {
            snpsWithin = new ArrayList<>(size());
            int len = getSequences().get(0).getLength();
            for (int k = 0; k < getSequences().size() - 1; k++) {
                MsaSequence s1 = getSequences().get(k);
                if (s1.getNonTipPaddingLength() > 0) {
                    continue;
                }
                boolean[] padding1 = s1.getPaddingArray(maxIndelLength, minIndelDistFromEnds);
                for (int l = 1; l < getSequences().size(); l++) {
                    MsaSequence s2 = getSequences().get(l);
                    if (s2.getNonTipPaddingLength() > 0) {
                        continue;
                    }
                    boolean[] padding2 = s2.getPaddingArray(maxIndelLength, minIndelDistFromEnds);
                    for (int i = 0; i < len; i++) {
                        if (s1.getSequenceString().charAt(i) != s2.getSequenceString().charAt(i) && !padding1[i] && !padding2[i]) {
                            snpsWithin.add(new Snp(s1, s2, i));
                        }
                    }
                }

            }
        }
    }
}
