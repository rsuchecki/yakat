/*
 * Copyright 2017 Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au.
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
package kextender;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.BlockingQueue;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class PairMersExtenderProducer implements Runnable {

    private final PairMersMap map;
    private final BlockingQueue<List<PairMer>> queue;
    private final int BUFFER_SIZE;

    public PairMersExtenderProducer(PairMersMap map, BlockingQueue<List<PairMer>> outqueue, int BUFFER_SIZE) {
        this.map = map;
        this.queue = outqueue;
        this.BUFFER_SIZE = BUFFER_SIZE;
    }

    @Override
    public void run() {

        try {
            Iterator<PairMer> it = map.getTerminalPairMers().keySet().iterator();
            ArrayList<PairMer> buffer = new ArrayList<>(BUFFER_SIZE);
            while (it.hasNext()) {
                PairMer pairMer = it.next();
                if (!pairMer.isVisited()) {
                    if (buffer.size() >= BUFFER_SIZE) {
                        putOneQueue(queue, buffer);
                        buffer = new ArrayList<>();
                    }
                    buffer.add(pairMer);
                }
            }
            if (!buffer.isEmpty()) {
                putOneQueue(queue, buffer);                
            }
//            for (int i = 0; i < consumers; i++) {
                putOneQueue(queue, new ArrayList<>(0)); //inform other threads                
//            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
//        pairMersMaps.addToTotalPairMersGenerated(pairMersGenerated);
    }

    private void putOneQueue(BlockingQueue q, List buffer) throws InterruptedException {
//        System.err.println(Thread.currentThread().getName() + "[P] puts " + buffer.size() + " on "+q.hashCode());
        q.put(buffer);
//        System.err.println(Thread.currentThread().getName() + "[P] puts " + buffer.size()+ " on "+q.hashCode() + " done");
    }
}
