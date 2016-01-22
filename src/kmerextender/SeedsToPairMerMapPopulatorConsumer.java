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

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SeedsToPairMerMapPopulatorConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> queue;
    private final PairMersMap pairMersMap;
    private final Integer KMER_LENGTH;

    public SeedsToPairMerMapPopulatorConsumer(BlockingQueue<ArrayList<String>> queue, PairMersMap pairMersMap,
            boolean splitInputSequenceintoKmers, Integer k, Integer minFreq) {
        this.queue = queue;
        this.pairMersMap = pairMersMap;
        this.KMER_LENGTH = k;
    }

    @Override
    public void run() {
        try {
            ArrayList<String> list;
            while (!(list = queue.take()).isEmpty()) {
                for (String read : list) {
                    ArrayList<String> kmers = getKmers(read.toUpperCase());
                    for (String kmerString : kmers) {
                        addKmerToMap(kmerString, KMER_LENGTH - 1);
                    }
                }
            }
            queue.put(new ArrayList<String>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private boolean addKmerToMap(String kmerString, int overlapLength) {
        try {
            pairMersMap.addToPairMersMap(kmerString, true, overlapLength, false);
            pairMersMap.addToPairMersMap(kmerString, false, overlapLength, false);
        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error !", this.getClass().getName());
            pairMersMap.setOutOfMemory();
            try {
                queue.put(new ArrayList<String>());
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
            return false;
        }
        return true;
    }

    /**
     * TODO    ---- ENCODE DIRECTLY TO PAIRMERS
     * @param sequence
     * @return 
     */
    private ArrayList<String> getKmers(String sequence) {
        int totalKmers = sequence.length() - KMER_LENGTH + 1;
        if (totalKmers < 0) { //in case of short sequences
            return new ArrayList<>(0);
        }
        ArrayList<String> kmers = new ArrayList<>(totalKmers);
        for (int i = 0; i < totalKmers; i++) {
            String s = sequence.substring(i, i + KMER_LENGTH);
            kmers.add(s);
        }
        return kmers;
    }

}
