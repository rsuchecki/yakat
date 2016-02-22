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
public class PairMersExtender {

    private final String DEBUG_FILE;
    private final String STATS_FILE;
    private final String TOOL_NAME;
//    private HashMap<String,String> ID_2_SEED_SEQUENCE_MAP;

    /**
     * Given optional DEBUG and STATS filenames, init Object
     *
     * @param DEBUG_FILE
     * @param STATS_FILE
     * @param TOOL_NAME
     */
    public PairMersExtender(String DEBUG_FILE, String STATS_FILE, String TOOL_NAME) {
        this.DEBUG_FILE = DEBUG_FILE;
        this.STATS_FILE = STATS_FILE;
        this.TOOL_NAME = TOOL_NAME;

    }

//    /**
//     * Given optional DEBUG and STATS filenames and optional seed(s), init Object
//     *
//     * @param DEBUG_FILE
//     * @param STATS_FILE
//     * @param idToSeedSequenceMap
//     */
//    public PairMersExtender(String DEBUG_FILE, String STATS_FILE, HashMap<String,String> idToSeedSequenceMap) {
//        this.DEBUG_FILE = DEBUG_FILE;
//        this.STATS_FILE = STATS_FILE;
//        this.ID_2_SEED_SEQUENCE_MAP = idToSeedSequenceMap;
//    }
    public void matchAndExtendKmers(int k, PairMersMap pairMersMap, boolean outputFasta, String namePrefix) {
        if (!namePrefix.isEmpty() && !namePrefix.endsWith("_")) {
            namePrefix += "_";
        }
        NavigableSet<PairMer> pairMers = pairMersMap.getPairMersSkipListMap().keySet();
        Iterator<PairMer> it = pairMers.iterator();
        long clusterNumber = 0; //Connected-component in the de-bruijn graph
        long extendedNumber = 0;
        int[] extendedLengths = null;
        int longest = 0;
        int MAX_LENGTH_STATS = 2000;
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", "Kmer extrending stats", TOOL_NAME), false);
            extendedLengths = new int[MAX_LENGTH_STATS];
        }
        while (it.hasNext()) {
            PairMer pairMer = it.next();
            if (!pairMer.isVisited()) {
                clusterNumber++;
                ConnectedPairMers connectedPairMers = new ConnectedPairMers();
                connectedPairMers.connectPairMers(pairMer, k, pairMersMap);
                try {
                    if (connectedPairMers.hasTerminalOrSingletonNode()) {
                        extendedNumber++;
//                        System.err.println(connectedPairMers.getKeys().size() + " pairMers == " + connectedMers.length() + " bp");
                        String connectedMers = connectedPairMers.toString(k);
                        int len = connectedMers.length();
                        if (len > longest) {
                            longest = len;
                        }
                        if (STATS_FILE != null) {
                            if (len < MAX_LENGTH_STATS) {
                                extendedLengths[len]++;
                            } else {
                                extendedLengths[0]++;
                            }
                        }
                        if (outputFasta) {
                            System.out.println(">" + namePrefix + clusterNumber + " " + len);
                            System.out.println(connectedMers);
                        } else {
                            System.out.println(connectedMers);// + "\t" + SequenceOps.getReverseComplementString(connected));
                        }

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
        String longestExtMessage = "Longest extended sequence " + NumberFormat.getNumberInstance().format(longest) + "bp";
        String totalExtendedMessage = "Total extended sequences = " + NumberFormat.getNumberInstance().format(extendedNumber);
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", longestExtMessage, TOOL_NAME), true);
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", totalExtendedMessage, TOOL_NAME), true);
            Reporter.writeHistogramToFile(STATS_FILE, extendedLengths, true, k);
        }
        Reporter.report("[INFO]", longestExtMessage, TOOL_NAME);
        Reporter.report("[INFO]", totalExtendedMessage, TOOL_NAME);
    }

//    public void matchAndExtendSeeds(int k, PairMersMap pairMersMap, PairMerToSeedMap pairMerToSeedMap) {
//        long clusterNumber = 0; //Connected-component in the de-bruijn graph
//
//        //Iterate through PairMers (2 per seed, corresponding to it's ends)
//        Iterator<PairMer> it = pairMerToSeedMap.getPairMersSkipListMap().keySet().iterator();
//        while (it.hasNext()) {
//            PairMer seedMer = it.next();
//            PairMer pairMer = pairMersMap.get(seedMer);
//            if (pairMer != null && !pairMer.isVisited()) {
//                clusterNumber++;
//                ConnectedPairMers connectedPairMers = new ConnectedPairMers();
//                connectedPairMers.connectPairMers(pairMer, k, pairMersMap);
//                try {
//                    if (connectedPairMers.hasTerminalOrSingletonNode()) {
//                        SeedSequence seedSequence = pairMerToSeedMap.get(seedMer);
//                        if (seedSequence != null && seedSequence.getSequenceString() != null) {
//                            String connectedMers = connectedPairMers.toString(k);
//                            String connectedMersRC = SequenceOps.getReverseComplementString(connectedMers);
////                            String extension = extendSeed(connectedMers, extendSeed(connectedMersRC, seedSequence.getExtendedOrOriginal(k), k), k);
////                            seedSequence.setExtended(k, extension);
//                            
//                            //experimenting with storing left and right extensions separately
//                            extendSeedAndSetExtensions(connectedMers, seedSequence, k);
//                            extendSeedAndSetExtensions(connectedMersRC, seedSequence, k);
//                        }
//                    } else {
//                        String message = "No terminal PairMer identified in cluster " + clusterNumber + " @ k=" + k;
//                        Reporter.report("[WARNING]", message, TOOL_NAME);
//                        if (DEBUG_FILE != null) {
//                            ArrayList<String> toReport = new ArrayList<>();
//                            toReport.add("No terminal PairMer identified in cluster " + clusterNumber + " @ k=" + k + " PairMers in cluster:");
//                            for (PairMer pm : connectedPairMers.getKeys()) {
//                                toReport.add(pm.getPairMerString(k));
//                            }
//                            Reporter.writeToFile(DEBUG_FILE, toReport, true);
//                        }
//                    }
//                } catch (StackOverflowError error) {
//                    Reporter.report("[ERROR]", "StackOverflow error, possible solution lies in : 'java -Xss<size> : set java thread stack size'", TOOL_NAME);
//                    Reporter.report("[ERROR]", "Obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", TOOL_NAME);
//                    System.exit(1);
//                }
//            }
//        }
//
//    }
//
//    private String extendSeed(String connectedMers, String seed, int k) {
//        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
//        String connectedMersHead = connectedMers.substring(0, k);
//        if (seed.endsWith(connectedMersHead)) {
////            System.err.println(seed);
////            System.err.println("Extended with ");
////            System.err.println(connectedMers);
//            return seed + connectedMers.substring(k);
//        }
//        String connectedMersTail = connectedMers.substring(connectedMers.length() - k, connectedMers.length());
//        if (seed.startsWith(connectedMersTail)) {
////            System.err.println(seed);
////            System.err.println("Extended with ");
////            System.err.println(connectedMers);
//            return connectedMers.substring(0, connectedMers.length() - k) + seed;
//        }
//        return seed;
//    }
//    
//    private void extendSeedAndSetExtensions(String connectedMers, SeedSequence seedSequence, int k) {
//        String seed = seedSequence.getSequenceString();
//        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
//        String connectedMersHead = connectedMers.substring(0, k);
//        if (seed.endsWith(connectedMersHead)) {
//            seedSequence.setRightExtension(k, connectedMers.substring(k));
//        }
//        String connectedMersTail = connectedMers.substring(connectedMers.length() - k, connectedMers.length());
//        if (seed.startsWith(connectedMersTail)) {
//            seedSequence.setLeftExtension(k, connectedMers.substring(0, connectedMers.length() - k));
//        }
//    }
//
////    private boolean checkIfExtensionPossible(ConnectedPairMers connectedPairMers, String seed, int k) {
////        //check both ends of the seed if k-1 end bases of the connected overlap either in forward or RC
//////        PairMer terminal1 = connectedPairMers.getTerminal1();
//////            for (PairMer terminal : new PairMer[]{connectedPairMers.getTerminal1(), connectedPairMers.getTerminal2()}) {
////        for (PairMer terminal : connectedPairMers.terminalMersOrSingleton()) {
////            String core = terminal.decodeCore(k - 1);
////            if (seed.endsWith(core) || seed.startsWith(core)) {
////                return true;
////            }
////            String rcCore = SequenceOps.getReverseComplementString(core);
////            if (seed.endsWith(rcCore) || seed.startsWith(rcCore)) {
////                return true;
////            }
////        }
////        return false;
//////        String connectedMersTail = connectedPairMers.substring(connectedPairMers.length() - k, connectedPairMers.length());
//////        if (seed.startsWith(connectedMersTail)) {
//////            return connectedPairMers.substring(0, connectedPairMers.length() - k) + seed;
//////        }
//////        return seed;
////    }
}
