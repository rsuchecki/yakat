/*
 * Copyright 2017 rad.
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
package snpmers;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListMap;
import shared.Sequence;
import shared.SequenceOps;

/**
 *
 * @author rad
 */
public class MapPopulatorConsumer implements Runnable {

    private final ConcurrentSkipListMap<String, ArrayList<KmerLink>> map;
    private final boolean DEBUG = false;
    private final BlockingQueue<ArrayList<SnpFilter>> inputQueue;
    private final int k;
    private final String TOOL_NAME;

    /**
     *
     * @param map
     * @param inputQueue
     * @param k
     * @param TOOL_NAME
     */
    public MapPopulatorConsumer(ConcurrentSkipListMap<String, ArrayList<KmerLink>> map, BlockingQueue<ArrayList<SnpFilter>> inputQueue,
        int k, String TOOL_NAME) {
        this.map = map;
        this.inputQueue = inputQueue;
        this.k = k;
        this.TOOL_NAME = TOOL_NAME;
    }

    @Override
    public void run() {

        try {
            ArrayList<SnpFilter> buffer = null;
            String previous = null;
            while (!(buffer = inputQueue.take()).isEmpty()) {
                for (SnpFilter snpFilter : buffer) {
                    kmerizeAndAddToMap(snpFilter, k, map);
                }

            }
//            inputQueue.put(new LabelledInputBuffer(null, new ArrayList())); //inform other threads
            inputQueue.put(new ArrayList()); //inform other threads
//            Reporter.report("[INFO]", "Finished assigning k-mer frequencies to SNPs", TOOL_NAME);

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void kmerizeAndAddToMap(SnpFilter snpFilter, int k, ConcurrentSkipListMap<String, ArrayList<KmerLink>> map) { //,
//            HashMap<String, ArrayList<KmerLink>> nonUniqueLinks) {
        //TWO PARENT SEQUENCES FOR EACH SNP
        for (int parent = 1; parent < 3; parent++) {
            Sequence parentSequence = snpFilter.getParentSequence(parent);
            String sequence = parentSequence.getSequenceString();
            int snpSite = snpFilter.getSnpPosition0();
            int offset = 0; //If padding shifts snpSite, store that offset here
            int maxKmer = Math.min(sequence.length() - k + 1, snpSite + 1);

            if (DEBUG) {
                for (int i = 0; i < sequence.length(); i++) {
                    if (i == snpSite) {
                        System.err.print("|");
                    } else {
                        System.err.print("_");
                    }
                }
                System.err.println();
                System.err.println(sequence);
            }

            //CANNOT K-MERIZE GAPS INSERTED BY MSA, SO NEED TO REMOVE THEM AND ADJUST THE SNP POSITION ACCORDINGLY
            if (sequence.contains("-")) {
                StringBuilder unpaddedSeq = new StringBuilder();
                int unpaddedSnpPosition = snpSite;
                for (int i = 0; i < sequence.length(); i++) {
                    if (sequence.charAt(i) == '-') {
                        if (i < snpSite) {
                            --unpaddedSnpPosition;
                            offset++;
                        }
                    } else {
                        unpaddedSeq.append(sequence.charAt(i));
                    }
                }
                sequence = unpaddedSeq.toString();
                snpSite = unpaddedSnpPosition;
                maxKmer = Math.min(sequence.length() - k + 1, snpSite + 1);

                if (DEBUG) {
                    for (int i = 0; i < sequence.length(); i++) {
                        if (i == snpSite) {
                            System.err.print("|");
                        } else {
                            System.err.print("_");
                        }
                    }
                    System.err.println();
                    System.err.println(sequence);
                }
            }

            int startAt = Math.max(0, snpSite - k + 1);
            for (int i = startAt; i < maxKmer; i++) {
                CharSequence kmer = sequence.subSequence(i, i + k);
                String canonical = SequenceOps.getCanonical(kmer.toString());

                int pos = i + offset; //position in the original/padded MSA sequence
//                KmerLink kmerLink = new KmerLink(snpFilter, (parent == 1), pos, !kmer.equals(canonical));
                KmerLink kmerLink = new KmerLink(snpFilter, (parent == 1), pos);
                ArrayList<KmerLink> kmerLinks = new ArrayList();
                kmerLinks.add(kmerLink);
                ArrayList<KmerLink> previousStored = map.putIfAbsent(canonical, kmerLinks);
                if (previousStored != null) {
                    synchronized (previousStored) {
                        previousStored.add(kmerLink);
                    }
                }                
                snpFilter.incrementMerCount(kmerLink.isParentOne());

//                ArrayList<KmerLink> kmerLinks = map.get(canonical);
//                if (kmerLinks == null) {
//                    kmerLinks = new ArrayList();
//                }
//                kmerLinks.add(kmerLink);
//                snpFilter.incrementMerCount(kmerLink.isParentOne());
//                map.put(canonical, kmerLinks);
//                if (map.containsKey(canonical)) {
//                    kmerLink.setUnique(false);
//                }
//                ArrayList<KmerLink> put = map.put(canonical, kmerLink);
//                if (put != null) {
//                    String key = parentSequence.getId() + "," + put.getParentSequence().getId();
//
//                    ArrayList<KmerLink> nonUniqList = nonUniqueLinks.getOrDefault(key, new ArrayList<KmerLink>());
//                    if (!nonUniqList.contains(kmerLink)) {
//                        nonUniqList.add(kmerLink);
//                    }
//                    if (!nonUniqList.contains(put)) {
//                        nonUniqList.add(put);
//                    }
//                    nonUniqueLinks.put(key, nonUniqList);
////                    Reporter.report("[WARNING]", "Non-unique ["+parentSequence.getId()+","+put.getParentSequence().getId()+"] k-mer overlapping a SNP will be ignored: " + canonical, TOOL_NAME);
////                        Reporter.report("[ERROR]", "Unexpected issue of a k-mer already present in the map: " + canonical, TOOL_NAME);
//                    if (DEBUG) {
//                        System.err.println("kmer, canonical:");
//                        System.err.println(kmer + ", " + canonical);
////                        System.err.println("Current:\n" + parentSequence.getId() + " at0=" + startPosition);
//                        System.err.println("Current:\n" + parentSequence.getId() + " at0=" + pos);
//                        System.err.println(parentSequence.getSequenceString());
//                        System.err.println("Previous:\n" + put.getParentSequence().getId() + " at0=" + put.getStartPosition());
//                        System.err.println(put.getParentSequence().getSequenceString());
//                    }
//                }
                if (DEBUG) {
                    for (int j = 0; j < i; j++) {
                        System.err.print(" ");
                    }
                    System.err.println(kmer);
                }
            }
            if (DEBUG) {
                System.err.println();
            }

        }
    }
}
