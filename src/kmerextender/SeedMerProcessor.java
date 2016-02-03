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

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SeedMerProcessor implements Runnable {

//    private final BlockingQueue<PairMersMap> queue;
    
    //Common 
    private final SeedSequences seedSequences;
    ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeeds;
    private final Integer MIN_KMER_FREQUENCY;
    private final String TOOL_NAME;
    
    
    private final PairMersMap pairMersMap;
    private final Integer KMER_SIZE;

    /**
     * 
     * @param seedSequences
     * @param kToSeeds
     * @param MIN_KMER_FREQUENCY
     * @param TOOL_NAME
     * @param pairMersMap
     * @param KMER_SIZE 
     */
    public SeedMerProcessor(SeedSequences seedSequences, ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeeds, Integer MIN_KMER_FREQUENCY, String TOOL_NAME, PairMersMap pairMersMap, Integer KMER_SIZE) {
        this.seedSequences = seedSequences;
        this.kToSeeds = kToSeeds;
        this.MIN_KMER_FREQUENCY = MIN_KMER_FREQUENCY;
        this.TOOL_NAME = TOOL_NAME;
        this.pairMersMap = pairMersMap;
        this.KMER_SIZE = KMER_SIZE;
    }

    

    @Override
    public void run() {
//        try {
//            PairMersMap map;
//            while (!(map = queue.take()).isEmpty()) {
//                long purged = map.purge();
//                if (purged > 100000) {
//                    gc(5, 500); //force GC 
//                }
//            }
//            queue.put(new PairMersMap()); //inform other threads

            PairMerToSeedMap pairMerToSeedMap = new PairMerToSeedMap(seedSequences, KMER_SIZE);
            PairMerToSeedMap previous = kToSeeds.putIfAbsent(KMER_SIZE, pairMerToSeedMap);
            if(previous != null) {
                Reporter.report("[BUG!]", "Overwritten exisitng value in "+kToSeeds.getClass().getSimpleName(), TOOL_NAME);
            }
            BlockingQueue<ArrayList<String>> seedsDummyQueue = new ArrayBlockingQueue<>(2);
            PairMersMap trimmedSeedPairMersMap = new PairMersMap();
            try {
                seedsDummyQueue.put(seedSequences.getSeedSequenceStrings());
                seedsDummyQueue.put(new ArrayList<String>());
            } catch (InterruptedException ex) {
            }
            PairMerMapPopulatorConsumer seedsPairMerMapPopulatorConsumer = new PairMerMapPopulatorConsumer(seedsDummyQueue,
                trimmedSeedPairMersMap, true, KMER_SIZE, MIN_KMER_FREQUENCY, true);
            seedsPairMerMapPopulatorConsumer.run(); 
            Reporter.report("[INFO]", "Seed-mers map populated, k=" + KMER_SIZE, TOOL_NAME);
            Iterator<PairMer> it = trimmedSeedPairMersMap.getPairMersSkipListMap().keySet().iterator();
            long removedCount = 0L;
            while (it.hasNext()) {
                if (pairMersMap.getPairMersSkipListMap().remove(it.next()) != null) {
                    removedCount++;
                }
            }
            Reporter.report("[INFO]", "Seed-based purging removed additional " + NumberFormat.getIntegerInstance().format(removedCount), TOOL_NAME);

//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }
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
