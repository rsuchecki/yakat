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
package kmerloc;

import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.atomic.AtomicLong;
import kmermatch.Kmer;
import kmermatch.KmerSetsMap;
import shared.Reporter;
import shared.Sequence;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class KmerLocatorConsumer implements Runnable {

    private final BlockingQueue<ArrayList<Sequence>> inputQueue;
    private final String TOOL_NAME;
    private final PrintStream bufferedOut;
    private final boolean storeASCII;
    private final KmerSetsMap kmerSetsMap;
//    private final LocatorStats stats;
//    private final BlockingQueue<Stat> statsQ;
    private final AtomicLong total;

    public KmerLocatorConsumer(BlockingQueue<ArrayList<Sequence>> inputQueue, KmerSetsMap kmerSetsMap, String TOOL_NAME, PrintStream bufferedOut,
            boolean storeASCII, AtomicLong total) {
//            boolean storeASCII, BlockingQueue<Stat> statsQ) {
//            boolean storeASCII, LocatorStats stats) {
        this.inputQueue = inputQueue;
        this.TOOL_NAME = TOOL_NAME;
        this.bufferedOut = bufferedOut;
        this.storeASCII = storeASCII;
        this.kmerSetsMap = kmerSetsMap;
        this.total = total;
//        this.stats = stats;
//        this.statsQ = statsQ;
//        this.k = k;
    }

    @Override
    public void run() {
        ArrayList<Sequence> list;
        try {

            while (!(list = inputQueue.take()).isEmpty()) {
                for (Sequence sequence : list) {
                    String id = sequence.getId();
                    ConcurrentHashMap<Integer, ConcurrentSkipListSet<Kmer>> map = kmerSetsMap.getKmerSetsMap();
                    Reporter.report("[INFO]", "Locating k-mers on " + sequence.getId(), TOOL_NAME);
                    long perSeqCount = 0;
                    for (Map.Entry<Integer, ConcurrentSkipListSet<Kmer>> entry : map.entrySet()) {
                        int count = 0;
                        int k = entry.getKey();
                        ConcurrentSkipListSet kmers = entry.getValue();
                        //KMERIZE
                        int maxKmer = sequence.getLength() - k + 1;
                        for (int i = 0; i < maxKmer; i++) {
                            CharSequence subSequence = sequence.getSequenceString().subSequence(i, i + k);                            
                            if (subSequence.chars().anyMatch(x -> x == 'N' || x == 'n')) {
                                continue;
                            }
                            String canonical = SequenceOps.getCanonical(subSequence.toString());
                            if (kmers.contains(new Kmer(canonical, storeASCII))) {
                                count++;
//                                String orientation = canonical.contentEquals(subSequence) ? "+" : "-";
//                                bufferedOut.append(id + ":" + (i + 1) + "-" + (i + k) + "\t" + canonical + System.lineSeparator());                                
                                bufferedOut.append(subSequence + "\t"+ id + "\t" + (i + 1) + "\t" + (i + k) +"\t"+System.lineSeparator());
                            }
                        }
//                        stats.addStat(id, k, count);
                        perSeqCount += count;                        
//                    Reporter.report("[INFO]", "Located " + count + " " + k + "-mers in " + id, TOOL_NAME);
                    }
                    Reporter.report("[INFO]", "Located " + NumberFormat.getInstance().format(perSeqCount) + " k-mers in " + id, TOOL_NAME);
                    total.getAndAdd(perSeqCount);
                }
            }
            inputQueue.put(new ArrayList<>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
