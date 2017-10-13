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
package snpmers;

import shared.BaseCall;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import shared.Reporter;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class SnpFilter implements Comparable<SnpFilter>{

    private final Sequence sequence1;
    private final Sequence sequence2;
    private final int snpPosition0;
//    private final int alignmentLength;
    private final String clusterId;
    private short[] mers1; //RECORD START POSITIONS OF ENCOUNTERED k-mers WHICH OVERLAP WITH THE SNP
    private short[] mers2; //RECORD START POSITIONS OF ENCOUNTERED k-mers WHICH OVERLAP WITH THE SNP
    private HashMap<String, BaseCall> snpCalls;
//    private HashMap<String, String> callDetails;
    private boolean valid = true;

    private int mersParent1;
    private int mersParent2;

    /**
     *
     * @param clusterId
     * @param sequence1
     * @param sequence2
     * @param snpPosition0 zero-indexed!
     * @param TOOL_NAME
     */
    public SnpFilter(String clusterId, Sequence sequence1, Sequence sequence2, int snpPosition0, String TOOL_NAME) {
        this.clusterId = clusterId;
        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
        this.snpPosition0 = snpPosition0;
        if(getBase1() == getBase2()) {
            Reporter.report("[ERROR]", "Input sequences in "+clusterId+" homozygous at reported SNP position: "+snpPosition0, TOOL_NAME);
            System.exit(1);
        }
        mers1 = new short[snpPosition0 + 1]; //
        mers2 = new short[snpPosition0 + 1];
    }

    public void clearMers(String sample) {
        mers1 = new short[snpPosition0 + 1]; //
        mers2 = new short[snpPosition0 + 1];
    }

    public String getClusterId() {
        return clusterId;
    }

    public boolean isIndel() {
        return getBase1() == '-' || getBase2() == '-';
    }

    public Sequence getParentSequence(int i) {
        if (i == 1) {
            return sequence1;
        } else if (i == 2) {
            return sequence2;
        }
        return null;
    }

    public Sequence getSequence1() {
        return sequence1;
    }

    public Sequence getSequence2() {
        return sequence2;
    }

    /**
     * Zero-indexed
     *
     * @return
     */
    public int getSnpPosition0() {
        return snpPosition0;
    }

    /**
     * Zero-indexed
     *
     * @return
     */
    public int getSnpPosition0UnpaddedSeq1() {
        String prefix = getSequence1().getSequenceString().substring(0, getSnpPosition0());
        String unpaddedPrefix = prefix.replaceAll("-", "");
        int diff = prefix.length() - unpaddedPrefix.length();
        return snpPosition0 - diff;
    }

    /**
     * Zero-indexed
     *
     * @return
     */
    public int getSnpPosition0UnpaddedSeq2() {
        String prefix = getSequence2().getSequenceString().substring(0, getSnpPosition0());
        String unpaddedPrefix = prefix.replaceAll("-", "");
        int diff = prefix.length() - unpaddedPrefix.length();
        return snpPosition0 - diff;
    }

    public boolean setMer1(int position, short value) {
        if (mers1[position] == 0) {
            mers1[position] = value;
            return true;
        }
        return false;
    }

    public boolean setMer2(int position, short value) {
        if (mers2[position] == 0) {
            mers2[position] = value;
            return true;
        }
        return false;
    }

    public short[] getMers1() {
        return mers1;
    }

    public short[] getMers2() {
        return mers2;
    }

    private ArrayList<Short> getNonZeroKmerFreqs(short[] mers) {
        ArrayList<Short> nonZeroValues = new ArrayList<>(mers.length);
        for (Short value : mers) {
            if (value > 0) {
                nonZeroValues.add(value);
            }
        }
        return nonZeroValues;
    }

    
    
    private double getMax(ArrayList<Short> nonZeroValues) {
        if (nonZeroValues.isEmpty()) {
            return 0;
        } else {
            short max = 0;
            for (Short nonZeroValue : nonZeroValues) {
                max = nonZeroValue > max ? nonZeroValue : max;
            }
            return max;
        }
    }

    public BaseCall getBaseCall(String id) {
        return snpCalls.get(id);
    }

//    public String getSnpCallDetails(String id) {
//        return callDetails.get(id);
//    }
    public void callBaseAndResetMers(String sampleName, int minTotal, int minMinor, double minCoverage, double maxError, String TOOL_NAME) {
        if (snpCalls == null) {
            snpCalls = new HashMap<>();
//            callDetails = new HashMap<>();
        }
//        System.out.println();
//        System.out.printf("%15s%30s%30s\n", clusterId, sequence1.getId(), sequence2.getId());
//        System.out.println(sequence1.getSequenceString());
//        System.out.println(sequence2.getSequenceString());
//        System.out.println("Calling "+sampleName+" parent calls "+getBase1()+"/"+getBase2());
        BaseCall call = call(minTotal, minMinor, minCoverage, maxError);
        BaseCall put = snpCalls.put(sampleName, call);

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
        if (put != null) {
            Reporter.report("[WARNING]", "Call " + put + " previously made for " + sampleName + ", current call: " + call, this.getClass().getSimpleName());
        }
//        if(this.clusterId.equals("Cluster_172")) {
//            int x = 0;
//        }
        mers1 = new short[mers1.length];
        mers2 = new short[mers2.length];
    }

//    private BaseCall call(int minTotal, int minMinor, int minKmers, String TOOL_NAME) {
    private BaseCall call(int minTotal, int minMinor, double minCoverageRatio, double maxError) {        
//        System.out.println(Arrays.toString(getMers1()).replaceAll(",|\\[|\\]", ""));
//        System.out.println(Arrays.toString(getMers2()).replaceAll(",|\\[|\\]", ""));
        
        ArrayList<Short> nonZeroKmerFreqs1 = getNonZeroKmerFreqs(getMers1());
        ArrayList<Short> nonZeroKmerFreqs2 = getNonZeroKmerFreqs(getMers2());
        double medianFreq1 = shared.CommonMaths.getMedian(nonZeroKmerFreqs1);
        double medianFreq2 = shared.CommonMaths.getMedian(nonZeroKmerFreqs2);
//        double cov1 = getMax(nonZeroKmerFreqs1);
//        double cov2 = getMax(nonZeroKmerFreqs2);

//        if(mersParent1 < 10 &&  mersParent2 < 10) {
//            int x =0;
//        }
        int cov1 = nonZeroKmerFreqs1.size();
        int cov2 = nonZeroKmerFreqs2.size();
        double coverageRatio1 = (double)cov1 / getMersParent1();
        double coverageRatio2 = (double)cov2 / getMersParent2();
        //QUICKLY DISCARD IF LOW FREQUENCY/LOW COVERAGE 
        if (medianFreq1 + medianFreq2 < minTotal || (coverageRatio1 < minCoverageRatio && coverageRatio2 < minCoverageRatio)) {
//            System.out.println(cov1+" "+cov2+" "+coverageRatio1+" "+coverageRatio2+" _/_ "+nonZeroKmerFreqs1.size()+"/"+getMersParent1()+" "+nonZeroKmerFreqs2.size()+"/"+getMersParent2());
            return new BaseCall(null, null, medianFreq1, medianFreq2, cov1, cov2);
        }
        Character base1 = getBase1();
        Character base2 = getBase2();
        //ALLOW FOR CERTAIN AMOUNT OF ERROR 
        //e.g. a few k-mers from elswhere in the genome could be ignored rather than leading to an ambigous call.
        if (coverageRatio1 <= maxError) {
            coverageRatio1 = 0;
            base1 = null;
        }
        if (coverageRatio2 <= maxError) {
            coverageRatio2 = 0;
            base2 = null;
        }
        //CLEAR-CUT HOMOZYGOUS CASES
        if (medianFreq1 == 0 && coverageRatio2 >= minCoverageRatio) { //homozygous
//        System.out.println(cov1+" "+cov2+" "+coverageRatio1+" "+coverageRatio2+" "+base2+"/_ "+nonZeroKmerFreqs1.size()+"/"+getMersParent1()+" "+nonZeroKmerFreqs2.size()+"/"+getMersParent2());
            return new BaseCall(base2, null, medianFreq1, medianFreq2, cov1, cov2);
        } else if (medianFreq2 == 0 && coverageRatio1 >= minCoverageRatio) { //homozygous
//            System.out.println(cov1+" "+cov2+" "+coverageRatio1+" "+coverageRatio2+" "+base1+"/_ "+nonZeroKmerFreqs1.size()+"/"+getMersParent1()+" "+nonZeroKmerFreqs2.size()+"/"+getMersParent2());
            return new BaseCall(base1, null, medianFreq1, medianFreq2, cov1, cov2);
        }
        
        if (medianFreq1 >= minMinor && medianFreq2 >= minMinor && coverageRatio1 >= minCoverageRatio && coverageRatio2 >= minCoverageRatio) {
//            System.out.println(cov1+" "+cov2+" "+coverageRatio1+" "+coverageRatio2+" "+base1+" "+base2 +" "+nonZeroKmerFreqs1.size()+"/"+getMersParent1()+" "+nonZeroKmerFreqs2.size()+"/"+getMersParent2());
            return new BaseCall(base1, base2, medianFreq1, medianFreq2, cov1, cov2);
        }
//        System.out.println(cov1+" "+cov2+" "+coverageRatio1+" "+coverageRatio2+" _/_ "+nonZeroKmerFreqs1.size()+"/"+getMersParent1()+" "+nonZeroKmerFreqs2.size()+"/"+getMersParent2());
        return new BaseCall(null, null, medianFreq1, medianFreq2, cov1, cov2);

    }

    public char getBase1() {
        return getSequence1().getSequenceString().charAt(getSnpPosition0());
    }

    public char getBase2() {
        return getSequence2().getSequenceString().charAt(getSnpPosition0());
    }

    public boolean isValid() {
        return valid;
    }

    public void setInvalid() {
        valid = false;
    }

    public void incrementMerCount(boolean parentOne) {
        if (parentOne) {
            mersParent1++;
        } else {
            mersParent2++;
        }
    }

    public int getMersParent1() {
        return mersParent1;
    }

    public int getMersParent2() {
        return mersParent2;
    }

    @Override
    public int compareTo(SnpFilter o) {
        int compareTo = getClusterId().compareTo(o.getClusterId());
        if(compareTo == 0) {
            return getSnpPosition0() - o.getSnpPosition0();
        }
        return compareTo;
    }
    
    @Override
    public boolean equals(Object o) {
        return compareTo((SnpFilter) o) == 0;
    }

    @Override
    public int hashCode() {
        return getClusterId().hashCode()+getSnpPosition0();
    }

    
}
