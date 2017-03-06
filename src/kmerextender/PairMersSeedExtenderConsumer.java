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

import shared.Reporter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.Sequence;
import shared.SequenceOps;

/**
 * Given a Map of purged PairMers, identify and extend connected k-mers
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersSeedExtenderConsumer implements Runnable {

    private final String DEBUG_FILE;
//    private final String STATS_FILE;
    private final String TOOL_NAME;
    private final BlockingQueue<PairMersMap> pairMersMapsQueue;
    private final ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers;
    private final byte threadId;
//    private HashMap<String,String> ID_2_SEED_SEQUENCE_MAP;

    public PairMersSeedExtenderConsumer(BlockingQueue<PairMersMap> pairMersMapsQueue,
        ConcurrentHashMap<Integer, PairMerToSeedMap> kToSeedMers, String DEBUG_FILE, String TOOL_NAME, byte threadId) {
        this.DEBUG_FILE = DEBUG_FILE;
//        this.STATS_FILE = STATS_FILE;
        this.TOOL_NAME = TOOL_NAME;
        this.pairMersMapsQueue = pairMersMapsQueue;
        this.kToSeedMers = kToSeedMers;
        this.threadId = threadId;
    }

//    public void matchAndExtendSeeds(int k, PairMersMap pairMersMap, PairMerToSeedMap pairMerToSeedMap) {
    @Override
    public void run() {
        PairMersMap pairMersMap;
        try {
            while (!(pairMersMap = pairMersMapsQueue.take()).isNull()) {
                int k = pairMersMap.getK();
                //Iterate through PairMers (2 per seed, corresponding to it's ends)
                Iterator<PairMer> it = kToSeedMers.get(k).getPairMersSkipListMap().keySet().iterator();
                while (it.hasNext()) {
                    PairMer seedMer = it.next();
                    PairMer pairMer = pairMersMap.get(seedMer);
                    if (pairMer != null && !pairMer.isVisited()) {
                        ConnectedPairMers connectedPairMers = new ConnectedPairMers();
                        if(!connectedPairMers.connectPairMers(pairMer, k, pairMersMap, threadId, DEBUG_FILE, false)) {
                            System.err.println("Failed connecting pairmers "+this.getClass().getCanonicalName());
                        }
                        try {
                            if (connectedPairMers.hasTerminalOrSingletonNode()) {
                                SeedSequence seedSequence = kToSeedMers.get(k).get(seedMer);
                                 if (seedSequence != null && seedSequence.getSequenceString() != null) {
                                    String connectedMers = connectedPairMers.toCharSeq(k).toString(); //TODO CONSIDER SWITCHING TO CHAr SEQ
                                    String connectedMersRC = SequenceOps.getReverseComplementString(connectedMers);
//                            String extension = extendSeed(connectedMers, extendSeed(connectedMersRC, seedSequence.getExtendedOrOriginal(k), k), k);
//                            seedSequence.setExtended(k, extension);

                                    //experimenting with storing left and right extensions separately
                                    extendSeedAndStoreExtensions(connectedMers, seedSequence, k);
                                    extendSeedAndStoreExtensions(connectedMersRC, seedSequence, k);
                                }
                            } else {
                                String message = "No terminal PairMer identified in cluster @ k=" + k;
                                Reporter.report("[WARNING]", message, TOOL_NAME);
                                if (DEBUG_FILE != null) {
                                    ArrayList<String> toReport = new ArrayList<>();
                                    toReport.add("No terminal PairMer identified in cluster @ k=" + k + " PairMers in cluster:");
                                    for (PairMer pm : connectedPairMers.getKeys()) {
                                        toReport.add(pm.getPairMerString(k));
                                    }
                                    Reporter.writeToFile(DEBUG_FILE, toReport, true);
                                }
                            }
                        } catch (StackOverflowError error) {
                            Reporter.report("[ERROR]", "StackOverflow error, possible solution lies in : 'java -Xss<size> : set java thread stack size'", TOOL_NAME);
                            Reporter.report("[ERROR]", "Obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", TOOL_NAME);
                            System.exit(1);
                        }
                    }
                }
                 Reporter.report("[INFO]", "Finished extending for k=" + k, TOOL_NAME);
            }
            pairMersMapsQueue.put(new PairMersMap(null));
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

    }

    private String extendSeed(String connectedMers, String seed, int k) {
        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
        String connectedMersHead = connectedMers.substring(0, k);
        if (seed.endsWith(connectedMersHead)) {
//            System.err.println(seed);
//            System.err.println("Extended with ");
//            System.err.println(connectedMers);
            return seed + connectedMers.substring(k);
        }
        String connectedMersTail = connectedMers.substring(connectedMers.length() - k, connectedMers.length());
        if (seed.startsWith(connectedMersTail)) {
//            System.err.println(seed);
//            System.err.println("Extended with ");
//            System.err.println(connectedMers);
            return connectedMers.substring(0, connectedMers.length() - k) + seed;
        }
        return seed;
    }

    private void extendSeedAndStoreExtensions(String connectedMers, SeedSequence seedSequence, int k) {
        String seed = seedSequence.getSequenceString();
        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
        String connectedMersHead = connectedMers.substring(0, k);
        if (seed.endsWith(connectedMersHead)) {
            seedSequence.setRightExtension(k, connectedMers.substring(k));
        }
        String connectedMersTail = connectedMers.substring(connectedMers.length() - k, connectedMers.length());
        if (seed.startsWith(connectedMersTail)) {
            seedSequence.setLeftExtension(k, connectedMers.substring(0, connectedMers.length() - k));
        }
    }

}
