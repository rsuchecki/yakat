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
package kmerextender;

import java.util.concurrent.ConcurrentSkipListMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class PairMerToSeedMap {

    private final ConcurrentSkipListMap<PairMer, SeedSequence> pairMerToSeedMap;
    private final String TOOL_NAME;

    public PairMerToSeedMap(SeedSequences seedSequences, int k, String TOOL_NAME) {
        this.TOOL_NAME = TOOL_NAME;
        pairMerToSeedMap = new ConcurrentSkipListMap<>();
        for (SeedSequence s : seedSequences.getSeedSequences()) {
            if (s.getSequenceString().length() >= k) {
                addToSeedPairMersMap(s.getSequenceString(), 0, k - 1, false, s);
//                addToSeedPairMersMap(s.getSequenceString().substring(0, k), false, k - 1, s);
                int len = s.getSequenceString().length();
                addToSeedPairMersMap(s.getSequenceString(), len - k, len - 1, true, s);
//                addToSeedPairMersMap(s.getSequenceString().substring(len - k, len), true, k - 1, s);
            }
        }
    }

    public boolean isEmpty() {
        return pairMerToSeedMap == null;
    }

    /**
     *
     * @param charSequence
     * @param kmerFrom
     * @param kmerTo inclusive
     * @param frontClip
     * @param seedSequence
     */
    public final void addToSeedPairMersMap(CharSequence charSequence, int kmerFrom, int kmerTo, boolean frontClip, SeedSequence seedSequence) {//, boolean inputKmersUnique) {

//        System.err.println("Generating: "+charSequence.subSequence(kmerFrom, kmerTo+1)+" from "+seedSequence.getSequenceString()); //.subSequence(kmerFrom, kmerTo));
//        PairMer pairMer = PairMerGenerator.generatePairMer(kmerString, frontClip, overlapLength);
        PairMer pairMer = null;
        boolean failed = false;
        try {
            pairMer = PairMerTypeSelector.generatePairMer(charSequence, kmerFrom, kmerTo, frontClip, (byte)1);
        } catch (NonACGTException ex) {
//            Reporter.report("[WARNING]", ex.getMessage(), getClass().getCanonicalName());
            failed = true;
        }
        if (!failed) {
            //Atomic operation START
            SeedSequence previousStored = pairMerToSeedMap.putIfAbsent(pairMer, seedSequence);
            //Atomic operation END
            int k = kmerTo - kmerFrom + 1;
//        System.err.println("Generated:  "+pairMer.getPairMerString(k, "_")+" at k="+k);
            if (previousStored != null) {
                //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
//            previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique);

                pairMerToSeedMap.put(pairMer, new SeedSequence());
                StringBuilder message = new StringBuilder("Removing non-unique PairMer2seed at k=");
                message.append(kmerTo - kmerFrom + 1).append(": ");
                String previousId;
                if ((previousId = previousStored.getId()) != null) {
                    message.append(previousId).append(", ");
                }
                message.append(seedSequence.getId());
                message.append(", ").append(pairMer.getPairMerString(kmerTo - kmerFrom + 1, "_"));
                Reporter.report("[WARNING]", message.toString(), TOOL_NAME);
            }
        }
    }
//    /**
//     * r
//     *
//     * @param kmerString
//     * @param frontClip
//     * @param overlapLength : normally k-1
//     * @param seedSequence
//     * @param inputKmersUnique : true when a list of unique k-mers given as
//     * input, false if k-mers extracted directly from FASTA/FASTQ
//     */
//    public final void addToSeedPairMersMap(String kmerString, boolean frontClip, int overlapLength, SeedSequence seedSequence) {//, boolean inputKmersUnique) {
//    
//    /**
//     * 
//     * @param kmerString
//     * @param kmerFrom
//     * @param kmerTo inclusive
//     * @param frontClip
//     * @param seedSequence 
//     */
////    public final void addToSeedPairMersMap(CharSequence kmerString, int kmerFrom, int kmerTo, boolean frontClip, SeedSequence seedSequence) {//, boolean inputKmersUnique) {
//
//        System.err.println("Generating: "+kmerString+" from "+seedSequence.getSequenceString()); //.subSequence(kmerFrom, kmerTo));
//        PairMer pairMer = PairMerGenerator.generatePairMer(kmerString, frontClip, overlapLength);
//        System.err.println("Generated:  "+pairMer.getPairMerString(kmerString.length(), "_")+" at k="+kmerString.length());
////        PairMer pairMer = PairMerGenerator.generatePairMer( kmerString, kmerFrom, kmerTo, frontClip);
//        //Atomic operation START
//        SeedSequence previousStored = pairMerToSeedMap.putIfAbsent(pairMer, seedSequence);
//        //Atomic operation END
//        if (previousStored != null) {
//            //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
////            previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique);
//            System.err.println("Should not have happened in NIKS pipeline" + getClass().getCanonicalName());
//        }
//    }

//    /**
//     * Given a core string, retrieves the matching PairMer
//     *
//     * @param core, which will be converted to its canonical form
//     * @param k
//     * @return PairMer if present in Map, null otherwise
//     */
//    public SeedSequence get(String core, int k) {
//        return pairMerToSeedMap.get(PairMerGenerator.getPairMer(core, k));
//    }
    /**
     * Retrieve a PairMer with a matching core
     *
     * @param anotherPairMer
     * @return
     */
    public SeedSequence get(PairMer anotherPairMer) {
        return pairMerToSeedMap.get(anotherPairMer);
    }

    /**
     *
     * @return the underlying data structure
     */
    public ConcurrentSkipListMap getPairMersSkipListMap() {
        return pairMerToSeedMap;
    }

    public int size() {
        return this.pairMerToSeedMap.size();
    }
}
