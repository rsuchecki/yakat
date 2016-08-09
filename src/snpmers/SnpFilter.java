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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import shared.Reporter;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class SnpFilter {

    private final Sequence sequence1;
    private final Sequence sequence2;
    private final int snpPosition;
    private final String clusterId;
    private short[] mers1; //RECORD START POSITIONS OF ENCOUNTERED k-mers WHICH OVERLAP WITH THE SNP
    private short[] mers2; //RECORD START POSITIONS OF ENCOUNTERED k-mers WHICH OVERLAP WITH THE SNP
    private HashMap<String, Character> snpCalls;
    private HashMap<String, String> callDetails;
    private boolean valid = true; 

    /**
     *
     * @param sequence1
     * @param sequence2
     * @param snpPosition zero-indexed!
     */
    public SnpFilter(String clusterId, Sequence sequence1, Sequence sequence2, int snpPosition) {
        this.clusterId = clusterId;
        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
        this.snpPosition = snpPosition;
        mers1 = new short[snpPosition + 1]; //
        mers2 = new short[snpPosition + 1];
    }

    public void clearMers(String sample) {
        mers1 = new short[snpPosition + 1]; //
        mers2 = new short[snpPosition + 1];
    }

    public String getClusterId() {
        return clusterId;
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
     * @return 
     */
    public int getSnpPosition() {
        return snpPosition;
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

    private double getMedian(ArrayList<Short> nonZeroValues) {
        if (nonZeroValues.isEmpty()) {
            return 0;
        } else {
            Collections.sort(nonZeroValues);
            int size = nonZeroValues.size();
            if (size == 1) {
                return (double) nonZeroValues.get(0);
            } else if (size % 2 == 0) {
                int middle = size / 2;
                int sum = nonZeroValues.get(middle - 1) + nonZeroValues.get(middle);
                return ((double) (sum)) / 2;
            } else {
                int middle = size / 2;
                return (double) nonZeroValues.get(middle);
            }
        }
    }


    public Character getSnpCall(String id) {
        return snpCalls.get(id);
    }

    public String getSnpCallDetails(String id) {
        return callDetails.get(id);
    }

    public char callBaseAndResetMers(String sampleName, int minTotal, int minMinor, int minKmers, String TOOL_NAME) {
        if (snpCalls == null) {
            snpCalls = new HashMap<>();
            callDetails = new HashMap<>();
        }
        char call = call(minTotal, minMinor, minKmers, TOOL_NAME);
        Character put = snpCalls.put(sampleName, call);

        //DEBUGGING ONLY
        ArrayList<Short> nonZeroKmerFreqs1 = getNonZeroKmerFreqs(getMers1());
        ArrayList<Short> nonZeroKmerFreqs2 = getNonZeroKmerFreqs(getMers2());
        double median1 = getMedian(nonZeroKmerFreqs1);
        double median2 = getMedian(nonZeroKmerFreqs2);
        if (median1 > 0 || median2 > 0) {
            StringBuilder callDetail = new StringBuilder("[");
            if (median1 > 0) {
                callDetail.append(getSequence1().getSequenceString().charAt(getSnpPosition())).append(":").append((int) Math.ceil(median1));
                callDetail.append("*").append(nonZeroKmerFreqs1.size());
            }
            if (median2 > 0) {
                if (median1 > 0) {
                    callDetail.append("/");
                }
                callDetail.append(getSequence2().getSequenceString().charAt(getSnpPosition())).append(":").append((int) Math.ceil(median2));
                callDetail.append("*").append(nonZeroKmerFreqs2.size());
            }
            callDetail.append("]");
            callDetails.put(sampleName, callDetail.toString());
        } else {
            callDetails.put(sampleName, "");
        }

        //DEBUGGING ONLY
        if (put != null) {
            Reporter.report("[WARNING]", "Call " + put + " previously made for " + sampleName + ", current call: " + call, this.getClass().getSimpleName());
        }
        if(this.sequence1.getId().equals("2545_156383")) {
            int x = 0;
        }
        mers1 = new short[mers1.length];
        mers2 = new short[mers2.length];
        return call;
    }

    private char call(int minTotal, int minMinor, int minKmers, String TOOL_NAME) {
        ArrayList<Short> nonZeroKmerFreqs1 = getNonZeroKmerFreqs(getMers1());
        ArrayList<Short> nonZeroKmerFreqs2 = getNonZeroKmerFreqs(getMers2());
        double median1 = getMedian(nonZeroKmerFreqs1);
        double median2 = getMedian(nonZeroKmerFreqs2);
//        if(median1 == 6 && median2 == 2) {
//            int x = 0;
//        }
        if (median1 + median2 < minTotal || (nonZeroKmerFreqs1.size() < minKmers && nonZeroKmerFreqs2.size() < minKmers)) {
            return 'N';
        } else if (median1 == 0 && nonZeroKmerFreqs2.size() >= minKmers) { //homozygous
            return getSequence2().getSequenceString().charAt(getSnpPosition());
        } else if (median2 == 0 && nonZeroKmerFreqs1.size() >= minKmers) { //homozygous
            return getSequence1().getSequenceString().charAt(getSnpPosition());
        } else {
            if (median1 >= minMinor && median2 >= minMinor && nonZeroKmerFreqs1.size() >= minKmers && nonZeroKmerFreqs2.size() >= minKmers) {
                char base1 = getBase1();
                char base2 = getBase2();
                if (base1 == '-' || base2 == '-') {
//                    Reporter.report("[WARNING]", "Unable to call IUPAC code for: " 
//                        + getSequence1().getId() + ":" + base1 + ", " + getSequence2().getId() + ":" + base2, TOOL_NAME);
                    return 'N';
                }
                return getIupacCode(base1, base2);
            }
            return 'N';
        }
    }

    public char getBase1() {
        return getSequence1().getSequenceString().charAt(getSnpPosition());
    }

    public char getBase2() {
        return getSequence2().getSequenceString().charAt(getSnpPosition());
    }

    private char getIupacCode(char c1, char c2) {
        if (c1 > c2) {
            char tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        if (c1 == 'A' && c2 == 'T') {
            return 'W';
        } else if (c1 == 'C' && c2 == 'G') {
            return 'S';
        } else if (c1 == 'A' && c2 == 'C') {
            return 'M';
        } else if (c1 == 'G' && c2 == 'T') {
            return 'K';
        } else if (c1 == 'A' && c2 == 'G') {
            return 'R';
        } else if (c1 == 'C' && c2 == 'T') {
            return 'Y';
        } else {
            Reporter.report("[FATAL] ", "Error calling IUPAC code for: " + c1 + "," + c2, this.getClass().getSimpleName());
            return '0';
        }
    }

    public boolean isValid() {
        return valid;
    }

    public void setInvalid() {
        valid = false;
    }
    
    
}
