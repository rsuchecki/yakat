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

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import shared.Reporter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
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
    private int WRITER_BUFFER_SIZE = 8192;

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
    public void matchAndExtendKmers(int k, PairMersMap pairMersMap, boolean outputFasta, String namePrefix, String outFile, int minLen) {
        if (!namePrefix.isEmpty() && !namePrefix.endsWith("_")) {
            namePrefix += "_";
        }

//        Iterator<PairMer> it = pairMersMap.getPairMersMap().keySet().iterator();
        Iterator<PairMer> it = pairMersMap.getTerminalPairMers().keySet().iterator();
        long clusterNumber = 0; //Connected-component in the de-bruijn graph
        long extendedNumber = 0;
        long extendedLength = 0;
        long longEnough = 0;
        long longEnoughBp = 0;
        int[] extendedLengths = null;
        int longest = 0;
        int shortest = Integer.MAX_VALUE;
        int MAX_LENGTH_STATS = 2000;
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", "Kmer extrending stats", TOOL_NAME), false);
            extendedLengths = new int[MAX_LENGTH_STATS];
        }
        try {
            BufferedWriter out;
            if (outFile.endsWith(".gz")) {
                out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile)), "UTF-8"), WRITER_BUFFER_SIZE);
            } else {
                out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"), WRITER_BUFFER_SIZE);
            }
            while (it.hasNext()) {
                PairMer pairMer = it.next();
                
                
                StringBuilder otherCore = new StringBuilder();
                 String decodedCore = pairMer.decodeCore(k - 1);
                 if(pairMer.hasBothClips()) {
                     Reporter.report("[WARNING]", "Both clips present in a suposedly terminal pairmer!", getClass().getCanonicalName());
                 } else if(pairMer.hasLeftClip()) {
                    otherCore.append(pairMer.getClipLeft());
                    otherCore.append(decodedCore.subSequence(0, decodedCore.length() - 1));                     
                 } else if(pairMer.hasRightClip()) {
                    otherCore.append(decodedCore.subSequence(1, decodedCore.length()));                     
                    otherCore.append(pairMer.getClipRight());                    
                 } 
           

                
                 
                try {
                    it.remove();
                    PairMer otherTerminal = pairMersMap.getTerminal(otherCore, k);
                    if(otherTerminal != null)  {
                        //skipping one-base extension
                        pairMersMap.removeTerminal(otherTerminal);
                        continue;
                    }
                    pairMer = pairMersMap.get(otherCore, k);
                } catch (NonACGTException ex) {
                    Reporter.report("[WARNING]", "Unexpected NonACGTException caught", getClass().getCanonicalName());
                }

                if (pairMer != null && !pairMer.isVisited()) {
                    clusterNumber++;
                    ConnectedPairMers connectedPairMers = new ConnectedPairMers(DEBUG_FILE);
                    connectedPairMers.connectPairMers(pairMer, k, pairMersMap);
////                    if(DEBUG_FILE != null &&  connectedPairMers.isBug()) {
//                    if(DEBUG_FILE != null) {
//                        ArrayList<String> toReport = new ArrayList<>();
//                        toReport.add("Cluster " + clusterNumber + " @ k=" + k + " PairMers in cluster:");
//                        for (PairMer pm : connectedPairMers.getKeys()) {
//                            toReport.add(pm.getPairMerString(k));
//                        }
//                        Reporter.writeToFile(DEBUG_FILE, toReport, true);
//                        
//                    }
                    try {
                        if (connectedPairMers.hasTerminalOrSingletonNode()) {
                            extendedNumber++;
//                        System.err.println(connectedPairMers.getKeys().size() + " pairMers == " + connectedMers.length() + " bp");
                            String connectedMers = connectedPairMers.toString(k);
                            int len = connectedMers.length();
                            extendedLength += len;
                            if (len > longest) {
                                longest = len;
                            }
                            if (len < shortest) {
                                shortest = len;
                            }
                            if (STATS_FILE != null) {
                                if (len < MAX_LENGTH_STATS) {
                                    extendedLengths[len]++;
                                } else {
                                    extendedLengths[0]++;
                                }
                            }
                            if (len >= minLen) {
                                longEnough++;
                                longEnoughBp += len;
                                if (outputFasta) {
                                    out.write(">" + namePrefix + clusterNumber + " " + len);
                                    out.newLine();
                                }
                                out.write(connectedMers);// + "\t" + SequenceOps.getReverseComplementString(connected));
                                out.newLine();
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
                out.flush();
                }
            }
        } catch (UnsupportedEncodingException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        } catch (IOException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        }
        String longestExtMessage = "Longest extended sequence = " + NumberFormat.getNumberInstance().format(longest) + "bp";
        String totalExtendedMessage = "Number of extended sequences = " + NumberFormat.getNumberInstance().format(extendedNumber);
        String totalExtendedMessage2 = "Total length of extended sequences = " + NumberFormat.getNumberInstance().format(extendedLength) + " bp";
        String longEnoughMessage = "Number of reported sequences " + NumberFormat.getNumberInstance().format(shortest == Integer.MAX_VALUE ? minLen : shortest) + "bp and longer = "
                + NumberFormat.getNumberInstance().format(longEnough);
        String longEnoughMessage2 = "Total length of reported sequences " + NumberFormat.getNumberInstance().format(shortest == Integer.MAX_VALUE ? minLen : shortest) + "bp and longer = "
                + NumberFormat.getNumberInstance().format(longEnoughBp) + " bp";
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", longestExtMessage, TOOL_NAME), true);
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", totalExtendedMessage, TOOL_NAME), true);
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", longEnoughMessage, TOOL_NAME), true);
            Reporter.writeHistogramToFile(STATS_FILE, extendedLengths, true, k);
        }
        Reporter.report("[INFO]", longestExtMessage, TOOL_NAME);
        Reporter.report("[INFO]", totalExtendedMessage, TOOL_NAME);
        Reporter.report("[INFO]", totalExtendedMessage2, TOOL_NAME);
        Reporter.report("[INFO]", longEnoughMessage, TOOL_NAME);
        Reporter.report("[INFO]", longEnoughMessage2, TOOL_NAME);
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
