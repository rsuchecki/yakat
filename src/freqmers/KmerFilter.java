/*
 * Copyright 2016 rad.
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
package freqmers;

import java.util.HashMap;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class KmerFilter implements Comparable<KmerFilter>{
    private final String id;
    private final Sequence sequence1;
    private short[] mers; //RECORD START POSITIONS OF ENCOUNTERED k-mers 
    private final HashMap<String, KmerFilterStats> samplesToStatsMap;
    private int mersCount;
    private int maxUniqMers;
    
    
    /**
     *
     * @param sequence1
     * @param k
     * @param TOOL_NAME
     */
    public KmerFilter(String id, Sequence sequence1, int k,  String TOOL_NAME) {
        this.sequence1 = sequence1;
        mers = new short[this.sequence1.getLength()-k+1]; //
        this.id = id;
        samplesToStatsMap = new HashMap<>();
    }

    public void clearMers(String sample) {
        mers = new short[getSequence().getLength() + 1]; //
        mersCount = 0;
    }

    public Sequence getSequence() {
        return sequence1;
    }

    public boolean setMer(int position, short value) {
        if (mers[position] == 0) {
            mers[position] = value;
            mersCount++;
            return true;
        }
        return false;
    }

    public int getMaxUniqMers() {
        return maxUniqMers;
    }

    public void setMaxUniqMers(int maxUniqMers) {
        this.maxUniqMers = maxUniqMers;
    }

    

    public short[] getMers1() {
        return mers;
    }

    

//    public String getSnpCallDetails(String id) {
//        return callDetails.get(id);
//    }
    public void collectStatsAndResetMers(String sampleName, String TOOL_NAME) {
//        System.out.println();
//        System.out.printf("%15s%30s%30s\n", clusterId, sequence1.getId(), sequence2.getId());
//        System.out.println(sequence1.getSequenceString());
//        System.out.println(sequence2.getSequenceString());
//        System.out.println("Calling "+sampleName+" parent calls "+getBase1()+"/"+getBase2());
        
        //        System.out.println(Arrays.toString(getMers1()).replaceAll(",|\\[|\\]", ""));
//        System.out.println(Arrays.toString(getMers2()).replaceAll(",|\\[|\\]", ""));
        
        

//        ArrayList<Short> nonZeroKmerFreqs1 = getNonZeroKmerFreqs(getMers1());
//        double medianFreq1 = shared.CommonMaths.getMedian(nonZeroKmerFreqs1);
//        double cov1 = getMax(nonZeroKmerFreqs1);
//        double cov2 = getMax(nonZeroKmerFreqs2);

//        if(mersParent1 < 10 &&  mersParent2 < 10) {
//            int x =0;
//        }
//        int cov1 = nonZeroKmerFreqs1.size();
//        double coverageRatio1 = (double)cov1 / getmaxMers();
        
        
//        BaseCall put = snpCalls.put(sampleName, call);

//        //DEBUGGING ONLY
//        ArrayList<Short> nonZeroKmerFreqs1 = getNonZeroKmerFreqs(getMers1());
//        ArrayList<Short> nonZeroKmerFreqs2 = getNonZeroKmerFreqs(getMers2());
//
////        if (clusterId.equals("Cluster_548")) {
//        int size1 = nonZeroKmerFreqs1.size();
//        int size2 = nonZeroKmerFreqs2.size();
////        if (uniqeMersParent1 != size1 && uniqeMersParent2 != size2 && (size1 > 0 || size2 >0)) {
////            System.err.print("Calling "+sampleName+" "+clusterId);
////            System.err.printf(" %.2f\t%.2f", (double) size1 / uniqeMersParent1, (double) size2 / uniqeMersParent2);
////            System.err.println(", uniqmers counts max, obs: " + uniqeMersParent1 + " / " + uniqeMersParent2 + ", " + size1 + " / " + size2);
//        System.out.println(uniqeMersParent1 + "\t" + uniqeMersParent2);
////                        
////            System.err.println(sequence1.getSequenceString());
////            System.err.println(Arrays.toString(mers1));            
//////            for(int i=0; i<mers1.length; i++) {
//////                if(mers1[i] > 0) {
//////                    System.err.println("["+i+"] "+sequence1.getSequenceString().substring(i, i+45));
//////                }
//////            }            
////            System.err.println(sequence2.getSequenceString());
////            System.err.println(Arrays.toString(mers2));
//////            for(int i=0; i<mers2.length; i++) {
//////                if(mers2[i] > 0) {
//////                    System.err.println("["+i+"] "+sequence2.getSequenceString().substring(i, i+45));
//////                }
//////            }      
//            int x=0;
////        }   
//        }
//        }

//        double median1 = getMedian(nonZeroKmerFreqs1);
//        double median2 = getMedian(nonZeroKmerFreqs2);
//        if (median1 > 0 || median2 > 0) {
//            StringBuilder callDetail = new StringBuilder("[");
//            if (median1 > 0) {
//                callDetail.append(getSequence1().getSequenceString().charAt(getSnpPosition0())).append(":").append((int) Math.ceil(median1));
//                callDetail.append("*").append(nonZeroKmerFreqs1.size());
//            }
//            if (median2 > 0) {
//                if (median1 > 0) {
//                    callDetail.append("/");
//                }
//                callDetail.append(getSequence2().getSequenceString().charAt(getSnpPosition0())).append(":").append((int) Math.ceil(median2));
//                callDetail.append("*").append(nonZeroKmerFreqs2.size());
//            }
//            callDetail.append("]");
//            callDetails.put(sampleName, callDetail.toString());
//        } else {
//            callDetails.put(sampleName, "");
//        }
        //DEBUGGING ONLY
//        if (put != null) {
//            Reporter.report("[WARNING]", "Call " + put + " previously made for " + sampleName + ", current call: " + call, this.getClass().getSimpleName());
//        }
//        if(this.clusterId.equals("Cluster_172")) {
//            int x = 0;
//        }

        samplesToStatsMap.put(sampleName, new KmerFilterStats(mers, mersCount, getMaxUniqMers()));
        
        mers = new short[mers.length];
        mersCount =0;
    }


    public String getId() {
        return id;
    }

    public KmerFilterStats getKmerFilterStats(String sampleId) {
        return samplesToStatsMap.get(sampleId);
    }
    

    @Override
    public int compareTo(KmerFilter o) {
        return getId().compareTo(o.getId());
    }
    
    @Override
    public boolean equals(Object o) {
        return compareTo((KmerFilter) o) == 0;
    }

    @Override
    public int hashCode() {
        return getId().hashCode();
    }

    
}
