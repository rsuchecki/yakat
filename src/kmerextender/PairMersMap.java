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

import java.util.Iterator;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Wrapper around a concurrent collection (ConcurrentSkipListMap, but other
 * options possible) of PairMer objects. After the Map is populated, use
 * .purge()
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersMap extends shared.MerMap {

    private final ConcurrentSkipListMap<PairMer, PairMer> pairMersSkipListMap;

    private final int MAX_1LONG_ENCODE = 32; //2bits per nucl, signed long so should be 31, but can use sign bit if lex ordering not needed, so 32 allowed 
    private final int MAX_2LONG_ENCODE = 64;
    private final int MAX_3LONG_ENCODE = 96;
    private final int MAX_4LONG_ENCODE = 128;
    private final int MAX_5LONG_ENCODE = 160;
    

//    private boolean OutOfMemory;

    /**
     * Instantiate the Map
     */
    public PairMersMap() {
        pairMersSkipListMap = new ConcurrentSkipListMap<>();
    }

    /**
     * First tries to atomically add a k-mer to the Map, if this fails,
     * synchronized method is used to update the previously stored PairMer
     *
     * @param kmerString
     * @param frontClip
     * @param overlapLength : normally k-1
     * @param inputKmersUnique : true when a list of unique k-mers given as
     * input, false if k-mers extracted directly from FASTA/FASTQ
     */
    public void addToPairMersMap(String kmerString, boolean frontClip, int overlapLength, boolean inputKmersUnique) {
        SplitMer splitMer = new SplitMer(kmerString, frontClip, overlapLength);
//        System.err.println("Adding\t"+kmerString+"\t"+frontClip+"\t"+splitMer.getLeftClip()+"_"+splitMer.getCore()+"_"+splitMer.getRightClip());

        PairMer pairMer = splitMer.generatepairMer(kmerString);

        

        //Atomic operation START
        PairMer previousStoredPairMer = pairMersSkipListMap.putIfAbsent(pairMer, pairMer);
        //Atomic operation END
        if (previousStoredPairMer != null) {
            //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
            previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique);
        }
    }

    /**
     * Given a core string, retrieves the matching PairMer
     *
     * @param core, which will be converted to its canonical form
     * @param k
     * @return PairMer if present in Map, null otherwise
     */
    public PairMer get(String core, int k) {

        if (k - 1 <= MAX_1LONG_ENCODE) {
            return pairMersSkipListMap.get(new PairMer1LongEncoded(core));
        } else if (k - 1 <= MAX_2LONG_ENCODE) {
            return pairMersSkipListMap.get(new PairMer2LongEncoded(core));
        } else if (k - 1 <= MAX_3LONG_ENCODE) {
            return pairMersSkipListMap.get(new PairMer3LongEncoded(core));
        } else if (k - 1 <= MAX_4LONG_ENCODE) {
            return pairMersSkipListMap.get(new PairMer4LongEncoded(core));
        } else if (k - 1 <= MAX_5LONG_ENCODE) {
            return pairMersSkipListMap.get(new PairMer5LongEncoded(core));
        } else {
            return pairMersSkipListMap.get(new PairMerIntArrEncoded(core));
        }
    }

    /**
     *
     * @return the underlying data structure
     */
    public ConcurrentSkipListMap getPairMersSkipListMap() {
        return pairMersSkipListMap;
    }

    /**
     * Removes from the Map each PairMer which (i) represents ambiguous
     * extension (>2 k-mers matching the core) or (ii) has no extensions (only
     * one k-mer matching the core) or (iii) represents 2 conflicting k-mers
     * (matching core, but both extending in the same direction) Run only when
     * the set/map is fully populated.
     *
     * @param k : k-mer length
     * @return the number of elements removed
     */
    public long purge(int k) {
        long count = 0L;
        Iterator<PairMer> it = pairMersSkipListMap.keySet().iterator();

        int total = 0;
        while (it.hasNext()) {
            PairMer next = it.next();
            if (next.isInvalid() || next.getStoredCount() != 2) {
//                System.err.println("Removing:");//\n\t"+next.getClipLeft()+"-"+next.getTmpCore()+"-"+next.getClipRight());
//                for (String s : next.getHistory()) {
//                    System.err.println("\t" + s);
//                }
                pairMersSkipListMap.remove(next);
                count++;
//                System.err.println(next.getPairMerString(k)+" <- PURGED");
            }
        }
        
        
        //EXPERIMENTS only
        it = pairMersSkipListMap.keySet().iterator();
        int buckets = pairMersSkipListMap.size();
        byte[] hashKeys = new byte[buckets];
        while (it.hasNext()) {
            PairMer next = it.next();
            int hashcode = next.hashCode() % buckets;
            if (hashcode < 0) {
                hashKeys[-hashcode]++;
            } else {
                hashKeys[hashcode]++;
            }
            total++;
        }
//        double mean = (double) total / buckets;
//        double deviation = 0;
//        int single = 0;
//        int multi = 0;
//        int empty = 0;
//        for (int i = 0; i < hashKeys.length; i++) {
////            System.err.printf("%8d %8d\n",i,hashKeys[i]);
//            deviation += Math.pow((hashKeys[i] - mean), 2);
//            if (hashKeys[i] < 1) {
//                empty++;
//            } else if (hashKeys[i] == 1) {
//                single++;
//            } else {
//                multi++;
//            }
//        }
//        System.err.println("Mean  = " + mean);
//        System.err.println("StDev = " + Math.sqrt(deviation / buckets));
//        System.err.println("Empty = " + empty + " " + ((double) empty / buckets) + "%");
//        System.err.println("Single= " + single + " " + ((double) single / buckets) + "%");
//        System.err.println("Multi = " + multi + " " + ((double) multi / buckets) + "%");
////        for (int i = hashKeys.length-10000; i < hashKeys.length; i++) {
////            System.err.printf("%8d %8d\n",i,hashKeys[i]);
////            
////        }
        return count;
    }

//    public boolean isOutOfMemory() {
//        return OutOfMemory;
//    }
//
//    public synchronized void setOutOfMemory() {
//        this.OutOfMemory = true;
//    }
    
//    public int purgeSeedMers(PairMersMap seedMersMap) {
//        Iterator<PairMer> it = seedMersMap.getPairMersSkipListMap().keySet().iterator();
//        
//    }

}
