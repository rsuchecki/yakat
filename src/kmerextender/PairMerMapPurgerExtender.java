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

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerMapPurgerExtender implements Runnable {

    private final BlockingQueue<PairMersMap> queue;
    private final Integer KMER_LENGTH;

    public PairMerMapPurgerExtender(BlockingQueue<PairMersMap> queue, Integer KMER_LENGTH) {
        this.queue = queue;
        this.KMER_LENGTH = KMER_LENGTH;
    }

    @Override
    public void run() {
        try {
            PairMersMap map;
            while (!(map = queue.take()).isEmpty()) {
            }
            queue.put(new PairMersMap()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

}
