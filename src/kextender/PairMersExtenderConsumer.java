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
package kextender;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersExtenderConsumer implements Runnable {

    private final PairMersMap map;
    private final BlockingQueue<List<PairMer>> inqueue;
    private final BlockingQueue<List<ConnectedPairMers>> outqueue;
    private String DEBUG_FILE = null;
    private final String TOOL_NAME;
    private int k;
    private final int OUT_BUFFER_SIZE;
    private final byte threadId;
    private int extensionsTerminated = 0;

    public PairMersExtenderConsumer(PairMersMap map, BlockingQueue<List<PairMer>> inqueue, BlockingQueue<List<ConnectedPairMers>> outqueue,
            int k, String DEBUG_FILE, byte threadId, String TOOL_NAME, int BUFFER_SIZE) {
        this.map = map;
        this.inqueue = inqueue;
        this.outqueue = outqueue;
        this.DEBUG_FILE = DEBUG_FILE;
        this.k = k;
        this.threadId = threadId;
        this.TOOL_NAME = TOOL_NAME;
        this.OUT_BUFFER_SIZE = BUFFER_SIZE;
    }

    @Override
    public void run() {
        try {
            List<PairMer> list;
            ArrayList<ConnectedPairMers> outBuffer = new ArrayList<>(OUT_BUFFER_SIZE);
            while (!(list = inqueue.take()).isEmpty()) {
//                System.err.println(Thread.currentThread().getName() + "[PC] taken " + list.size());
                for (PairMer pairMer : list) {
                    ///extend, put on outqueue
//                                System.err.println(pairMer.getPairMerString(6, "_")+" AMB="+pairMer.isAmbiguousTMP());

                    if (!pairMer.isVisited()) {
                        try {
                            ConnectedPairMers connectedPairMers = new ConnectedPairMers();
                            if (!connectedPairMers.connectPairMers(pairMer, k, map, threadId, DEBUG_FILE)) {
                                extensionsTerminated++;
//                            Reporter.report("[INFO]", "Thread "+threadId+" duplicate extension detected for "+pairMer.getPairMerString(k, "_"), TOOL_NAME);
//                            System.err.println("Failed extending pairmer "+pairMer.getPairMerString(k));
                                continue;
                            }
//                        if(connectedPairMers.size() == 0) {
//                            int x = 0;
//                            connectedPairMers.connectPairMers(pairMer, k, map, threadId, DEBUG_FILE, true);
//                        }
                            connectedPairMers.toCharSeq(k);
                            
                            
                            
                            if (outBuffer.size() >= OUT_BUFFER_SIZE) {
                                putOneQueue(outqueue, outBuffer);
                                outBuffer = new ArrayList<>(OUT_BUFFER_SIZE);
                            }
                            outBuffer.add(connectedPairMers);
                        } catch (StackOverflowError error) {
                            Reporter.report("[ERROR]", "StackOverflow error, possible solution lies in : 'java -Xss<size> : set java thread stack size'", TOOL_NAME);
                            Reporter.report("[ERROR]", "To obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", TOOL_NAME);
                            System.exit(1);
                        }
                    }
                }
            }
            if (!outBuffer.isEmpty()) {
                putOneQueue(outqueue, outBuffer);
            }
            putOneQueue(outqueue, new ArrayList<>(0)); //inform other threads
            putOneQueue(inqueue, new ArrayList<>()); //inform other threads
//            Reporter.report("[INFO]", "Thread " + threadId + " finished, terminated extensions = " + extensionsTerminated, TOOL_NAME + " " + Thread.currentThread().getName());

        } catch (InterruptedException e) {
            e.printStackTrace();
        }
//        pairMersMaps.addToTotalPairMersGenerated(pairMersGenerated);
    }

    private void putOneQueue(BlockingQueue q, List buffer) throws InterruptedException {
//        System.err.println(Thread.currentThread().getName() + "[PC] puts " + buffer.size() + " on " + q.hashCode());
        q.put(buffer);
//        System.err.println(Thread.currentThread().getName() + "[PC] puts " + buffer.size() + " on " + q.hashCode() + " done");
    }

}
