/*
 * Copyright 2015 Australian Centre For Plant Functional Genomics
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

import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerToSeedMapPopulator implements Runnable {

//    private final BlockingQueue<PairMersMap> queue;
    //Common 
    private final SeedSequences seedSequences;
    ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeeds;
    private final String TOOL_NAME;

    private final Integer KMER_SIZE;

    /**
     *
     * @param seedSequences
     * @param kToSeeds
     * @param TOOL_NAME
     * @param k
     */
    public PairMerToSeedMapPopulator(SeedSequences seedSequences, ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeeds,
        String TOOL_NAME, Integer k) {
        this.seedSequences = seedSequences;
        this.kToSeeds = kToSeeds;
        this.TOOL_NAME = TOOL_NAME;
        this.KMER_SIZE = k;
    }

    @Override
    public void run() {
        PairMerToSeedMap pairMerToSeedMap = new PairMerToSeedMap(seedSequences, KMER_SIZE);
        PairMerToSeedMap previous = kToSeeds.putIfAbsent(KMER_SIZE, pairMerToSeedMap);
        if (previous != null) {
            Reporter.report("[BUG!]", "Overwritten exisitng value in " + kToSeeds.getClass().getSimpleName(), TOOL_NAME);
        }

    }

    private void gc(int iterations, int sleep) {
        for (int i = 0; i < iterations; i++) {
            System.gc();
            try {
                Thread.sleep(sleep);
            } catch (InterruptedException ex) {
            }
        }
    }
}
