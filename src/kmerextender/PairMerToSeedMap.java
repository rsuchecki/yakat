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

/**
 *
 * @author rad
 */
public class PairMerToSeedMap {

    private final ConcurrentSkipListMap<PairMer, SeedSequence> pairMerToSeedMap;

    public PairMerToSeedMap(SeedSequences seedSequences, int k) {
        pairMerToSeedMap = new ConcurrentSkipListMap<>();
        for (SeedSequence s : seedSequences.getSeedSequences()) {
            addToPairMersMap(s.getSequenceString().substring(0, k), false, k - 1, s);
            int len = s.getSequenceString().length();
            addToPairMersMap(s.getSequenceString().substring(len - k, len), true, k - 1, s);
        }
    }

    public boolean isEmpty() {
        return pairMerToSeedMap == null;
    }

    /**
     * r
     *
     * @param kmerString
     * @param frontClip
     * @param overlapLength : normally k-1
     * @param seedSequence
     * @param inputKmersUnique : true when a list of unique k-mers given as
     * input, false if k-mers extracted directly from FASTA/FASTQ
     */
    public final void addToPairMersMap(String kmerString, boolean frontClip, int overlapLength, SeedSequence seedSequence) {//, boolean inputKmersUnique) {

        PairMer pairMer = PairMerGenerator.generatePairMer(kmerString, frontClip, overlapLength);

        //Atomic operation START
        SeedSequence previousStored = pairMerToSeedMap.putIfAbsent(pairMer, seedSequence);
        //Atomic operation END
        if (previousStored != null) {
            //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
//            previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique);
            System.err.println("Should not have happened in NIKS pipeline" + getClass().getCanonicalName());
        }
    }

    /**
     * Given a core string, retrieves the matching PairMer
     *
     * @param core, which will be converted to its canonical form
     * @param k
     * @return PairMer if present in Map, null otherwise
     */
    public SeedSequence get(String core, int k) {
        return pairMerToSeedMap.get(PairMerGenerator.getPairMer(core, k));
    }

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
}
