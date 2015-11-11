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

/**
 * Given a Map of purged PairMers, identify and extend connected k-mers
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersExtender {

    private String DEBUG_FILE;
    private String STATS_FILE;

    /**
     * Given optional DEBUG and STATS filenames, init Object
     *
     * @param DEBUG_FILE
     * @param STATS_FILE
     */
    public PairMersExtender(String DEBUG_FILE, String STATS_FILE) {
        this.DEBUG_FILE = DEBUG_FILE;
        this.STATS_FILE = STATS_FILE;
    }

    public void matchAndExtendKmers(int k, PairMersMap pairMersMap, boolean outputFasta, String namePrefix, int threads) {
        NavigableSet<PairMer> pairMers = pairMersMap.getPairMersSkipListMap().keySet();

        Iterator<PairMer> it = pairMers.iterator();
        long clusterNumber = 0;
        long extendedNumber = 0;
        int[] extendedLengths = null;
        int longest = 0;
        int MAX_LENGTH_STATS = 2000;
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", "Kmer extrending stats", getClass().getSimpleName()), false);
            extendedLengths = new int[MAX_LENGTH_STATS];
        }
        long parallelizableTime=0;
        while (it.hasNext()) {
            PairMer pairMer = it.next();
//            System.err.println(pairMer.getClipLeft()+"_"+pairMer.getTmpCore()+"_"+pairMer.getClipRight()+" Current PairMer, visited="+pairMer.isVisited());
            if (!pairMer.isVisited()) {
                clusterNumber++;
//            ArrayList<String> connectedKmersStrings = new ArrayList<>();
//            PairMerWrapper merWrapper = new PairMerWrapper(pairMer, k);
//            connectedKmersStrings.add(merWrapper.getKmer1String());
//            connectedKmersStrings.add(merWrapper.getKmer2String());
                ConnectedPairMers connectedPairMers = new ConnectedPairMers();

                connectedPairMers.connectPairMers(pairMer, k, pairMersMap);
//                PairMer terminal = addPairMerToListOfConnected(pairMer, k, pairMersMap, connectedPairMers);
//                Reporter.report("[INFO]", connectedPairMers.size() + " connected pairMers in cluster " + clusterNumber);
//                if (terminal == null) {
//                    Reporter.report("[WARNING]", "No terminal PairMer identified in cluster " + clusterNumber);
                if (!connectedPairMers.hasTerminalOrSingletonNode()) {
                    String message = "No terminal PairMer identified in cluster " + clusterNumber;
                    Reporter.report("[WARNING]", message, getClass().getSimpleName());
                    if (DEBUG_FILE != null) {
                        Set<PairMer> keys = connectedPairMers.getKeys();
                        ArrayList<String> toReport = new ArrayList<>();
                        toReport.add("No terminal PairMer identified in cluster " + clusterNumber + " PairMers in cluster:");
                        for (PairMer pm : keys) {
                            toReport.add(pm.getPairMerString(k));
                        }
                        Reporter.writeToFile(DEBUG_FILE, toReport, true);
                    }
                } else {
                    extendedNumber++;
//                    Reporter.report("[INFO]", "Terminal: " + terminal.getPairMerString(k));
                    long time = System.nanoTime();
                    String connected = null;
                    try{
                        connected = connectedPairMers.toString(k);
                    } catch (StackOverflowError error) {
                        Reporter.report("[ERROR]", "StackOverflow error, possible solution lies in : 'java -Xss<size> : set java thread stack size'", getClass().getSimpleName());
                        Reporter.report("[ERROR]", "Obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", getClass().getSimpleName());
                        System.exit(1);
                    }
                    parallelizableTime += (System.nanoTime() - time);
                    int len = connected.length();
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
                        System.out.println(">" + namePrefix +"_" + clusterNumber+" "+len);
                        System.out.println(connected);
                    } else {
                        System.out.println(connected);// + "\t" + SequenceOps.getReverseComplementString(connected));
                    }
                }
//                System.err.println("Connected PairMers set: ");
//                for (PairMer pm : connectedPairMers.getKeys()) {
//                    System.err.println(pm.getPairMerString(k));
//                }
            }
        }
        String longestExtMessage = "Longest extended sequence "+NumberFormat.getNumberInstance().format(longest)+"bp";
        String totalExtendedMessage = "Total extended sequences = " + NumberFormat.getNumberInstance().format(extendedNumber);
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", longestExtMessage, getClass().getSimpleName()), true);
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", totalExtendedMessage, getClass().getSimpleName()), true);
            Reporter.writeHistogramToFile(STATS_FILE, extendedLengths, true, k);
        }        
        Reporter.report("[INFO]",  longestExtMessage, getClass().getSimpleName());
        Reporter.report("[INFO]", totalExtendedMessage, getClass().getSimpleName());
        Reporter.report("[INFO]", "Potentially parallelizable single-threaded time elapsed: "+ NumberFormat.getNumberInstance().format((double)parallelizableTime * 1E-9) + " seconds", getClass().getSimpleName());

    }

}
