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
package kmerextender;

import java.text.NumberFormat;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;
import static kmerextender.CoreCoder.decodeCore;
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class NewPairMerArrays {

    private PairMer[][] pairMerArrays;
    
    private final int prefixLen;
    private final int postfixLen;
    private String TOOL_NAME;
//    private BlockingQueue perPrefixQueue[];

    public NewPairMerArrays(int prefixLen, int postfixLen, String TOOL_NAME) {
        this.prefixLen = prefixLen;
        this.postfixLen = postfixLen;
        int prefixMax = (int) Math.pow(4, prefixLen);
        int postfixMax = (int) Math.pow(4, postfixLen);
        this.pairMerArrays = new PairMer[prefixMax][postfixMax];
        this.TOOL_NAME = TOOL_NAME;
//        this.perPrefixQueue = new BlockingQueue[prefixMax];
    }

    

    public boolean addMer(String core, char clip, boolean left, AtomicLongArray stats) throws InterruptedException {
        int prefix = CoreCoder.encodeCore(core.subSequence(0, prefixLen));
        int postfix = CoreCoder.encodeCore(core.subSequence(prefixLen, prefixLen + postfixLen));
        PairMer arr[] = pairMerArrays[prefix];
        NewPairMerIntArrEncoded newPairMerIntArrEncoded = null;
        try {
            newPairMerIntArrEncoded = new NewPairMerIntArrEncoded(core.subSequence(prefixLen + postfixLen, core.length()), clip, left, 1);
        } catch (NonACGTException ex) {
            ex.printStackTrace();
        }
//        perPrefixQueue[prefix].put(new NewPairMerWithPrefixes(prefix, postfix, newPairMerIntArrEncoded));
        synchronized (arr) {
            if (arr[postfix] == null) {
//                    arr[prefix][prefixp] = new PairMerIntArrEncoded(inputKmerSequence, prefixLen + postfixLen, inputKmerSequence.length()-1, false, 1);                                    
                arr[postfix] = newPairMerIntArrEncoded;
                stats.incrementAndGet(0);
                return true;
            } else {
                PairMer previous = arr[postfix];
                while (previous != null) {
//                    if (previous.decodeCore(core.length()).equals(core) && ((NewPairMerIntArrEncoded)previous).compareTo(newPairMerIntArrEncoded) != 0) {
//
//                        System.err.println(core + " supposedly diff from ");
//                        System.err.println(previous.decodeCore(core.length()) + " this one ");
//                        int x = 9;
//                    }
                    if (((NewPairMerIntArrEncoded) previous).compareTo(newPairMerIntArrEncoded) == 0) {
                        //same core - TODO - add PM to core                    
                        stats.incrementAndGet(1);
//                    stats[1]++;
                        return true;
                    } else {
//                        System.err.println(previous.decodeCore(coreReminder.length()) + " != " + coreReminder);
                    }
                    previous = previous.getNextPairMer();
                }
                //replace first elem in linked list
//                PairMer stored = arr[prefix][prefixp];                
//                NewPairMerIntArrEncodedWithRef toStore = new NewPairMerIntArrEncodedWithRef(coreReminder, clip, left, 1, arr[prefix][prefixp]);
                arr[postfix] = new NewPairMerIntArrEncodedWithRef(newPairMerIntArrEncoded.getBitFields(), clip, left, 1, arr[postfix]);
                stats.incrementAndGet(2);
                return true;
//                stats[2]++;
            }
        }
    }

//    public BlockingQueue[] getPerPrefixQueue() {
//        return perPrefixQueue;
//    }

    
    
    public void printStats() {
        //get some stats
        int lev1load = 0;
        long occupiedSlots = 0;
        long totalElements = 0;
        long singlePMs = 0;
        
        for (int i = 0; i < pairMerArrays.length; i++) {
            CharSequence prefixSeq = decodeCore(prefixLen, i);
            if (pairMerArrays[i] != null) {
                lev1load++;
//                int localCount = 0;
                for (int j = 0; j < pairMerArrays[i].length; j++) {
//                    CharSequence prefixpSeq = decodeCore(postfixLen, j);
                    if (pairMerArrays[i][j] != null) {
                        occupiedSlots++;
//                        localCount++;
                        int local = 0;                        
                        PairMer p = pairMerArrays[i][j];
                        do {
                            local++;
                        } while ((p = p.getNextPairMer()) != null);
//                        System.err.println(prefixSeq+"_"+prefixpSeq+"_"+arr[i][j].decodeCore(16-prefixLen-postfixLen)+" "+local);                        
                        totalElements += local;
                        if(local == 1 ) {
                            singlePMs++;
                        }
                    }
                }
//                System.err.println(prefixSeq + " total " + localCount);
            }
        }
        NumberFormat perc = NumberFormat.getPercentInstance();
        perc.setMaximumFractionDigits(4);
        NumberFormat numInstance = NumberFormat.getInstance();
        numInstance.setMaximumFractionDigits(2);
        int slotsAvailable = (int) Math.pow(4, prefixLen)*(int) Math.pow(4, postfixLen);
        Reporter.report("[INFO]", "Occupied level 1 = " + perc.format(((double) lev1load) / pairMerArrays.length), TOOL_NAME);
//        Reporter.report("[INFO]", "Occupied level 2 = " + perc.format(((double) lev2load / pairMerArrays[0].length / pairMerArrays.length)), TOOL_NAME);
        Reporter.report("[INFO]", "Mean non empty occupancy level 2 = " + numInstance.format(((double) totalElements / occupiedSlots)), TOOL_NAME);
        Reporter.report("[INFO]", "Slots available = " + numInstance.format(slotsAvailable), TOOL_NAME);
        Reporter.report("[INFO]", "Slots occupied = " + numInstance.format((occupiedSlots)) + " ("+perc.format((double)occupiedSlots/slotsAvailable)+")", TOOL_NAME);
        Reporter.report("[INFO]", "Slots occupied by single PM  = " + numInstance.format((singlePMs)) + " ("+perc.format((double)singlePMs/slotsAvailable)+")", TOOL_NAME);
        Reporter.report("[INFO]", "Total elements  = " + numInstance.format((totalElements)), TOOL_NAME);
    }

}
