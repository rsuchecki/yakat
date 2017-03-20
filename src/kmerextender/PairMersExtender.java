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
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import shared.InputReaderProducer;
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
    private final int WRITER_BUFFER_SIZE = 8192;

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
    public void matchAndExtendKmers(int k, PairMersMap pairMersMap, boolean outputFasta, String namePrefix, String outFile, int minLen, int extenderThreads) {
        if (!namePrefix.isEmpty() && !namePrefix.endsWith("_")) {
            namePrefix += "_";
        }

        long clusterNumber = 0; //Connected-component in the de-bruijn graph, not the output cluster - change that?
        long extendedNumber = 0;
        long extendedLength = 0;
        long longEnough = 0;
        long longEnoughBp = 0;
        int[] extendedLengths = null;
        int longest = 0;
        int shortest = Integer.MAX_VALUE;
        int MAX_LENGTH_STATS = 2000;
        int potentialDuplicate = 0;
        if (STATS_FILE != null) {
            Reporter.writeToFile(STATS_FILE, Reporter.formatReport("[STATS]", "Kmer extrending stats", TOOL_NAME), false);
            extendedLengths = new int[MAX_LENGTH_STATS];
        }
        BufferedWriter out = null;
        try {
            if (outFile.endsWith(".gz")) {
                out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile)), "UTF-8"), WRITER_BUFFER_SIZE);
            } else {
                out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"), WRITER_BUFFER_SIZE);
            }

            byte threadId = Byte.MIN_VALUE;
            BlockingQueue inqueue = new ArrayBlockingQueue(extenderThreads);
            BlockingQueue<List<ConnectedPairMers>> outqueue = new ArrayBlockingQueue(extenderThreads);
            ArrayList<Future<?>> futures = new ArrayList<>(extenderThreads + 1);
            final ExecutorService producerConsumerExecutor = new ThreadPoolExecutor(extenderThreads + 1, extenderThreads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
            //SPAWN PRODUCER
            int BUFFER_SIZE = Math.max(1, (int) Math.min(1024, pairMersMap.sizeTerminals()/extenderThreads/4));
            futures.add(producerConsumerExecutor.submit(new PairMersExtenderProducer(pairMersMap, inqueue, BUFFER_SIZE)));
            //SPAWN CONSUMER THREADS
            for (int i = 0; i < extenderThreads; i++) {
                PairMersExtenderConsumer consumer = new PairMersExtenderConsumer(pairMersMap, inqueue, outqueue, k, DEBUG_FILE, threadId++, TOOL_NAME);
                futures.add(producerConsumerExecutor.submit(consumer));
            }
            //COUNT AND PRINT EXTENDED SEQUENCES, RESOLVE CONFLICTS...
            try {
                List<ConnectedPairMers> list;
                long pickedFromOutqueue = 0;
                long reportThreshold = (long) 1000;

                while (!(list = outqueue.take()).isEmpty() || --extenderThreads > 0) {
//                    System.err.println(Thread.currentThread().getName() + "[M] taken " + list.size());
                    if(list.isEmpty()) { //A thread finished 
                        Reporter.report("[INFO]", "Extending thread finished...???", TOOL_NAME);
                        continue;
                    }
                    pickedFromOutqueue += list.size();
                    if (pickedFromOutqueue % reportThreshold == 0 && pickedFromOutqueue != 0) {
                        if (reportThreshold < 1e6) {
//                                    reportThreshold *= 10;
//                                    reportThreshold *= KMER_REPORTING_MULTIPLY;
                            reportThreshold <<= 1; // *= 2
                        }
                        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(pickedFromOutqueue) + " extended so far...", TOOL_NAME);
                    }
                    for (ConnectedPairMers connectedPairMers : list) {
                        clusterNumber++;
                        try {
                            if (connectedPairMers.hasTerminalOrSingletonNode()) {
                                extendedNumber++;
//                        System.err.println(connectedPairMers.getKeys().size() + " pairMers == " + connectedMers.length() + " bp");
                                CharSequence connectedMers = connectedPairMers.toCharSeq(k);
                                int len = connectedMers.length();
                                extendedLength += len;

                                if (connectedPairMers.isPotentialDuplicate()) {
                                    potentialDuplicate++;
                                }
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

                                    //SYNC THIS 
                                    if (outputFasta) {
                                        out.write(">" + namePrefix + clusterNumber + " " + len + (connectedPairMers.isPotentialDuplicate() ? " potential_duplicate" : ""));
                                        out.newLine();
                                    }
                                    out.write(connectedMers.toString());// + "\t" + SequenceOps.getReverseComplementString(connected));
                                    out.newLine();
//                            } else {
//                                System.err.println("too short at "+len+" "+connectedMers.toString()+" given minlen = "+minLen);
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
                            Reporter.report("[ERROR]", "To obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", TOOL_NAME);
                            System.exit(1);
                        }
                        out.flush();

                    }
//                    System.err.println(list.size()+" "+extenderThreads);
                }
                Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(pickedFromOutqueue) + " extended multithreaded", TOOL_NAME);
//                System.err.println("out of lock");
//                outqueue.put(new ArrayList<>()); //inform other threads

            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            producerConsumerExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                producerConsumerExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!", TOOL_NAME);
            } catch (ExecutionException ex) {
                Reporter.report("[ERROR]", "PairMerSet populator execution exception!", TOOL_NAME);
                ex.printStackTrace();

            } catch (TimeoutException ex) {
                Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
                Reporter.report("[ERROR]", "PairMerSet populator timeout exception!", TOOL_NAME);
            }

            Reporter.report("[INFO]", "Starting second-pass, single threaded extending", TOOL_NAME);
            Iterator<PairMer> it = pairMersMap.getPairMersMap().keySet().iterator();
            int extendedInSecondPass = 0;
            while (it.hasNext()) {
                PairMer pairMer = it.next();
                if (pairMer != null && !pairMer.isVisited()) {
                    clusterNumber++;
                    ConnectedPairMers connectedPairMers = new ConnectedPairMers();
                    if (!connectedPairMers.connectPairMers(pairMer, k, pairMersMap, threadId, DEBUG_FILE)) {
                        System.err.println("Failed connecting pairmers - shoukld not happen in second pass");
                    }
                    try {
                        if (connectedPairMers.hasTerminalOrSingletonNode()) {
                            extendedNumber++;
//                        System.err.println(connectedPairMers.getKeys().size() + " pairMers == " + connectedMers.length() + " bp");
                            CharSequence connectedMers = connectedPairMers.toCharSeq(k);
                            int len = connectedMers.length();
                            extendedLength += len;

                            if (len > longest) {
                                longest = len;
                            }
                            if (len < shortest) {
                                shortest = len;
                            }

//                            if (len != k + 1) {
//                                System.err.println(connectedMers);
//                                for (PairMer pm : connectedPairMers.getKeys()) {
//                                    System.err.println(pm.getPairMerString(k, "_"));
//                                }
//                            }
//                            
                            if (STATS_FILE != null) {
                                if (len < MAX_LENGTH_STATS) {
                                    extendedLengths[len]++;
                                } else {
                                    extendedLengths[0]++;
                                }
                            }
                            if (len >= minLen) {
                                extendedInSecondPass++;
                                longEnough++;
                                longEnoughBp += len;
                                if (outputFasta) {
                                    out.write(">" + namePrefix + clusterNumber + " " + len);
                                    out.newLine();
                                }
                                out.write(connectedMers.toString());// + "\t" + SequenceOps.getReverseComplementString(connected));
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
                        Reporter.report("[ERROR]", "To obtain defaults: 'java -XX:+PrintFlagsFinal -version | grep ThreadStackSize'", TOOL_NAME);
                        System.exit(1);
                    }
                }
            }
            out.close();
            Reporter.report("[INFO]", "Extended in second-pass = " + NumberFormat.getNumberInstance().format(extendedInSecondPass) + " (single threaded extending)", TOOL_NAME);

        } catch (UnsupportedEncodingException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        } catch (IOException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
//        } catch (InterruptedException e) {
//            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        } finally {
            try {
                if (out != null) {
                    out.close();
                }
            } catch (IOException ex) {
                Reporter.report("[ERROR]", ex.getMessage(), TOOL_NAME);                
            }

        }
        String longestExtMessage = "Longest extended sequence = " + NumberFormat.getNumberInstance().format(longest) + " bp";
        String totalExtendedMessage = "Number of extended sequences = " + NumberFormat.getNumberInstance().format(extendedNumber);
        String totalExtendedMessage2 = "Total length of extended sequences = " + NumberFormat.getNumberInstance().format(extendedLength) + " bp";
        String longEnoughMessage = "Number of reported sequences " + NumberFormat.getNumberInstance().format(minLen) + " bp or longer = "
                + NumberFormat.getNumberInstance().format(longEnough);
        String longEnoughMessage2 = "Total length of reported sequences " + NumberFormat.getNumberInstance().format(minLen) + " bp or longer = "
                + NumberFormat.getNumberInstance().format(longEnoughBp) + " bp";
        String duplicatesExtMessage = "Number of terminated extensions = " + NumberFormat.getNumberInstance().format(potentialDuplicate);
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
        Reporter.report("[WARNING]", duplicatesExtMessage, TOOL_NAME);

    }

}
