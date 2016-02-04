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

import java.util.concurrent.BlockingQueue;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerMapPurger implements Runnable {

    private final BlockingQueue<PairMersMap> queue;

    public PairMerMapPurger(BlockingQueue<PairMersMap> queue) {
        this.queue = queue;
    }

    @Override
    public void run() {
        try {
            PairMersMap map;
            while (!(map = queue.take()).isEmpty()) {
                long purged = map.purge();
                if (purged > 100000) {
                    gc(5, 500); //force GC 
                }
            }
            queue.put(new PairMersMap(null)); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
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
