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
import shared.Sequence;
import shared.SequenceOps;

/**
 * Given a Map of purged PairMers, identify and extend connected k-mers
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersSeedExtenderConsumer  implements Runnable {

    private final String DEBUG_FILE;
    private final String STATS_FILE;
    private final String TOOL_NAME;
    private final int k;
    private final boolean FASTA_OUT;
    private String NAME_PREFIX;
    private final PairMersMap pairMersMap;
    private final PairMerToSeedMap pairMerToSeedMap;
//    private HashMap<String,String> ID_2_SEED_SEQUENCE_MAP;

//    /**
//     * Given optional DEBUG and STATS filenames, init Object
//     *
//     * @param DEBUG_FILE
//     * @param STATS_FILE
//     * @param TOOL_NAME
//     */
//    public PairMersSeedExtenderConsumer(String DEBUG_FILE, String STATS_FILE, String TOOL_NAME,
//        int k, PairMersMap pairMersMap, boolean outputFasta, String namePrefix, PairMerToSeedMap pairMerToSeedMap) {
//        this.DEBUG_FILE = DEBUG_FILE;
//        this.STATS_FILE = STATS_FILE;
//        this.TOOL_NAME = TOOL_NAME;
//        
//    }

    public PairMersSeedExtenderConsumer(String DEBUG_FILE, String STATS_FILE, String TOOL_NAME, 
        int k, boolean FASTA_OUT, String NAME_PREFIX, PairMersMap pairMersMap, PairMerToSeedMap pairMerToSeedMap) {
        this.DEBUG_FILE = DEBUG_FILE;
        this.STATS_FILE = STATS_FILE;
        this.TOOL_NAME = TOOL_NAME;
        this.k = k;
        this.FASTA_OUT = FASTA_OUT;
        this.NAME_PREFIX = NAME_PREFIX;
        this.pairMersMap = pairMersMap;
        this.pairMerToSeedMap = pairMerToSeedMap;
    }

    

    @Override
    public void run() {
        if (!NAME_PREFIX.isEmpty() && !NAME_PREFIX.endsWith("_")) {
            NAME_PREFIX += "_";
        }
        long clusterNumber = 0; //Connected-component in the de-bruijn graph
        Iterator<PairMer> it = pairMerToSeedMap.getPairMersSkipListMap().keySet().iterator();
        while (it.hasNext()) {
            PairMer seedMer = it.next();
            PairMer pairMer = pairMersMap.get(seedMer);
            if (pairMer != null && !pairMer.isVisited()) {
                clusterNumber++;
                ConnectedPairMers connectedPairMers = new ConnectedPairMers();
                connectedPairMers.connectPairMers(pairMer, k, pairMersMap);
                try {
                    if (connectedPairMers.hasTerminalOrSingletonNode()) {
//                        ArrayList<PairMer> terminalMersOrSingleton = connectedPairMers.terminalMersOrSingleton();
//                        for (PairMer endMer : terminalMersOrSingleton) {
                            SeedSequence seedSequence = pairMerToSeedMap.get(seedMer);
                            if (seedSequence != null) {
                                String connectedMers = connectedPairMers.toString(k);
                                String extension = extendSeed(connectedMers, extendSeed(SequenceOps.getReverseComplementString(connectedMers), seedSequence.getExtendedOrOriginal(k), k), k);
                                seedSequence.setExtended(k, extension);
                            }
//                        }
                    } else {
                        String message = "No terminal PairMer identified in cluster " + clusterNumber + " @ k=" + k;
                        Reporter.report("[WARNING]", message, TOOL_NAME);
                        if (DEBUG_FILE != null) {
                            ArrayList<String> toReport = new ArrayList<>();
                            toReport.add("No terminal PairMer identified in cluster " + clusterNumber + " @ k=" + k + " PairMers in cluster:");
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

//    private boolean checkIfExtensionPossible(ConnectedPairMers connectedPairMers, String seed, int k) {
//        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
////        PairMer terminal1 = connectedPairMers.getTerminal1();
////            for (PairMer terminal : new PairMer[]{connectedPairMers.getTerminal1(), connectedPairMers.getTerminal2()}) {
//        for (PairMer terminal : connectedPairMers.terminalMersOrSingleton()) {
//            String core = terminal.decodeCore(k - 1);
//            if (seed.endsWith(core) || seed.startsWith(core)) {
//                return true;
//            }
//            String rcCore = SequenceOps.getReverseComplementString(core);
//            if (seed.endsWith(rcCore) || seed.startsWith(rcCore)) {
//                return true;
//            }
//        }
//        return false;
////        String connectedMersTail = connectedPairMers.substring(connectedPairMers.length() - k, connectedPairMers.length());
////        if (seed.startsWith(connectedMersTail)) {
////            return connectedPairMers.substring(0, connectedPairMers.length() - k) + seed;
////        }
////        return seed;
//    }
}
