/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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
package kmermatch;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListSet;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerSetPopulatorConsumer implements Runnable {

    private final ConcurrentSkipListSet<Kmer> kmers;
    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final Integer k;

    /**
     *
     * @param kmers
     * @param inputQueue
     * @param k - set to null if no need to kmerize input
     */
    public KmerSetPopulatorConsumer(ConcurrentSkipListSet<Kmer> kmers, BlockingQueue<ArrayList<String>> inputQueue, Integer k) {
        this.kmers = kmers;
        this.inputQueue = inputQueue;
        this.k = k;
    }

    @Override
    public void run() {
        try {
            List<String> list;
            if (k ==null) { //KMERS INPUT, NO NEED TO KMERIZE
                while (!(list = inputQueue.take()).isEmpty()) {
                    for (String line : list) {
                        kmers.add(new KmerASCII(SequenceOps.getCanonical(line)));
                    }
                }
            } else {  //KMERIZE
                while (!(list = inputQueue.take()).isEmpty()) {
                    for (String line : list) {
                        int maxKmer = line.length() - k + 1;
                        for (int i = 0; i < maxKmer; i++) {
                            kmers.add(new KmerASCII(SequenceOps.getCanonical(line.subSequence(i, i + k).toString())));
                        }
                    }
                }
            }
            inputQueue.put(new ArrayList<>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

}
