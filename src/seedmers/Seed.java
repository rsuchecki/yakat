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
package seedmers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import shared.Reporter;
import shared.BaseCall;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class Seed {

    private final String id;
    private final CharSequence sequence;
    private final int snpPosition;
    private short mers[][];
    private HashMap<String, BaseCall> snpCalls;

    public Seed(String id, CharSequence sequence, int k) {
        this.id = id;
        this.sequence = sequence;
        this.snpPosition = k - 1;
        mers = new short[4][snpPosition + 1];
    }

    public String getId() {
        return id;
    }

    public CharSequence getSequence() {
        return sequence;
    }

    public int getSnpPosition() {
        return snpPosition;
    }

    public boolean setMer(AltSeedLink seedLink, short value) {
//        char allele = seedLink.getBaseChar();
        short mrs[] = mers[seedLink.getBase()];
//        switch (allele) {
//            case 'A':
//                mrs = mers[0];
//                break;
//            case 'C':
//                mrs = mers[1];
//                break;
//            case 'G':
//                mrs = mers[2];
//                break;
//            case 'T':
//                mrs = mers[3];
//                break;
//            default: {
//                Reporter.report("[WARNING]", "Unrecognized based in kmer " + kmer + " - failed setting mer count...", this.getClass().getSimpleName());
//                return false;
//            }
//        }
        if (mrs[seedLink.getPosition()] == 0) {
            mrs[seedLink.getPosition()] = value;
            return true;
        }
        return false;
    }

//    public void callBaseAndResetMers(String sampleName, int minTotal, int minMinor, double minCoverage, double maxError, String TOOL_NAME) {
    public void callBaseAndResetMers(String sampleName, String TOOL_NAME) {
        if (snpCalls == null) {
            snpCalls = new HashMap<>();
//            callDetails = new HashMap<>();
        }

//        BaseCall call = call(minTotal, minMinor, minCoverage, maxError);
        BaseCall call = call();
        BaseCall put = snpCalls.put(sampleName, call);
        //DEBUGGING ONLY
        if (put != null) {
            Reporter.report("[WARNING]", "Call " + put + " previously made for " + sampleName + ", current call: " + call, this.getClass().getSimpleName());
        }
//        if(this.clusterId.equals("Cluster_172")) {
//            int x = 0;
//        }        
        mers = new short[4][snpPosition + 1];

    }

    public BaseCall getBaseCall(String id) {
        return snpCalls.get(id);
    }

//    private BaseCall call(int minTotal, int minMinor, double minCoverageRatio, double maxError) {        
    private BaseCall call() {
//        int minTotal = 10;
//        int minMinor = 2; //double minCoverageRatio, double maxError
////        ArrayList<Short> nonZeroKmerFreqsA = getNonZeroKmerFreqs(mersA);
////        ArrayList<Short> nonZeroKmerFreqsC = getNonZeroKmerFreqs(mersC);
////        ArrayList<Short> nonZeroKmerFreqsG = getNonZeroKmerFreqs(mersG);
////        ArrayList<Short> nonZeroKmerFreqsT = getNonZeroKmerFreqs(mersT);
//        double covA = shared.CommonMaths.getMedian(nonZeroKmerFreqsA);
//        double covC = shared.CommonMaths.getMedian(nonZeroKmerFreqsC);
//        double covG = shared.CommonMaths.getMedian(nonZeroKmerFreqsG);
//        double covT = shared.CommonMaths.getMedian(nonZeroKmerFreqsT);
////        double cov1 = getMax(nonZeroKmerFreqs1);
////        double cov2 = getMax(nonZeroKmerFreqs2);
////        int uniqeMersParent1 = getUniqeMersParent1();
////        int uniqeMersParent2 = getUniqeMersParent2();
////        if(mersParent1 < 10 &&  mersParent2 < 10) {
////            int x =0;
////        }
//        double coverageRatioA = (double) nonZeroKmerFreqsA.size() / mersA.length;
//        double coverageRatioC = (double) nonZeroKmerFreqsC.size() / mersC.length;
//        double coverageRatioG = (double) nonZeroKmerFreqsG.size() / mersG.length;
//        double coverageRatioT = (double) nonZeroKmerFreqsT.size() / mersT.length;

       
        
        int numKmers[] = new int[mers.length];
        double coverageMedian[] = new double[mers.length];
//        double coveredRatio[] = new double[mers.length];
        for (int i = 0; i < mers.length; i++) {
            short[] m = mers[i];
            ArrayList<Short> nonZeroKmerFreqs = getNonZeroKmerFreqs(m);
            numKmers[i] = nonZeroKmerFreqs.size();
            coverageMedian[i] = shared.CommonMaths.getMedian(nonZeroKmerFreqs);
//            coveredRatio[i] = (double) nonZeroKmerFreqs.size() / m.length;
//            double cov = shared.CommonMaths.getMedian(nonZeroKmerFreqs);
//            double covRatio = (double) nonZeroKmerFreqs.size() / m.length;
//            sb.append("\t");               
//            sb.append(nonZeroKmerFreqs.size());
//            sb.append("\t");               
//            sb.append(cov);               
//            sb.append("\t");               
//            sb.append(CommonMaths.round(covRatio,2));               
        }
        StringBuilder sb = new StringBuilder(id);        
        sb.append("\t").append(getSequence().charAt(getSnpPosition()));
        sb.append("\t").append(Arrays.toString(numKmers));
        sb.append("\t").append(Arrays.toString(coverageMedian));
//        sb.append("\t").append(Arrays.toString(coveredRatio));
        sb.append("\t").append(getIUPAC(numKmers));
        sb.append("\t").append(sequence);        
        System.out.println(sb.toString());

        //QUICKLY DISCARD IF LOW FREQUENCY/LOW COVERAGE 
//        if (cov1 + cov2 < minTotal || (coverageRatio1 < minCoverageRatio && coverageRatio2 < minCoverageRatio)) {
//            return new BaseCall(null, null);
//        }
//        Character base1 = getBase1();
//        Character base2 = getBase2();
//        //ALLOW FOR CERTAIN AMOUNT OF ERROR 
//        //e.g. a few k-mers from elswhere in the genome could be ignored rather than leading to an ambigous call.
//        if (coverageRatio1 <= maxError) {
//            coverageRatio1 = 0;
//            base1 = null;
//        }
//        if (coverageRatio2 <= maxError) {
//            coverageRatio2 = 0;
//            base2 = null;
//        }
//        //CLEAR-CUT HOMOZYGOUS CASES
//        if (cov1 == 0 && coverageRatio2 >= minCoverageRatio) { //homozygous
//            return new BaseCall(base2, null);
//        } else if (cov2 == 0 && coverageRatio1 >= minCoverageRatio) { //homozygous
//            return new BaseCall(base1, null);
//        }
//        
//        if (cov1 >= minMinor && cov2 >= minMinor && coverageRatio1 >= minCoverageRatio && coverageRatio2 >= minCoverageRatio) {
//            return new BaseCall(base1, base2);
//        }
        return new BaseCall(getIUPAC(numKmers));

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

    /**
     * Given an array of numbers representing observed base counts
     * [#A,#C,#G,#T]  call IUPAC representation of
     * the observed bases. Bases are not considered for IUPAC reporting if there
     * are less than minCoverageThreshold or
     *
     * @param encoded [#A,#C,#G,#T]
     * @param locusDepth
     * @param maxErrorPercent % of totalDepth threshold above which a base (or
     * resulting IUPAC) is reported rather than being ignored as erroneous
     * @paramminCoveragePerAlleleminCoverageThreshold int min depth for a base
     * to be considered
     *
     * TODO EXTEND TO ACCOMMODATE IDELS
     *
     * @return
     */
    private char getIUPAC(int[] encoded) {
//        if (locusDepth == 0) {
//            return '0';
//        } else if (locusDepth < minCoveragePerLocus || locusDepth > maxCoveragePerLocus) {
//            return '?';
//        }
        boolean tt[] = new boolean[4]; //truth table i=0,1,2,3 values A,C,G,T...
        //Convert max perc error into max int coverage
//        int maxErrAlleleInt = (int) Math.floor(locusDepth * maxPercErrAllele / 100);
//        System.err.println("MaxErrInt="+maxErrInt);
        //Max of maxErr min coverage threshold 
//        int minCoverage = Math.max((int) Math.ceil(totalDepth * (float) maxErrorPercent / 100), minCoverageThreshold);
        int total = 0;
        for (int i = 0; i < encoded.length; i++) {
            if (encoded[i] > 0) {
                tt[i] = true;
            }
//            if (encoded[i] < minCoverageThreshold || encoded[i] <= maxErrInt) {
//                totalOther += encoded[i];                
//            }

//            if (encoded[i] >= minCoveragePerAllele && encoded[i] > maxErrAlleleInt) {
//                tt[i] = true;
//            } else {
//                totalOther += encoded[i];
                total += encoded[i];
//            }
        }
        char base = '?'; //ambiguousCallChar;
//        if (totalOther > maxErrAlleleInt) {
//            return ambiguousCallChar;
//        }
        if (!tt[0] && !tt[1] && !tt[2] && !tt[3]) {
//            if (tt[5] && !tt[6]) { //insertion in reference
//                base = 'E';
//            } else if (!tt[5] && tt[6]) { //deletion in reference
//                base = 'I';
//            } else if (tt[5] && tt[6]) { //deletion and insertion??
//                base = ambiguousCallChar;
//            } else 
            if (total == 0) {
                base = '.';
//            } else {
//                base = '?';
            }
//        } else if ((tt[0] || tt[1] || tt[2] || tt[3]) && (tt[5] || tt[6])) { //base and an indel???
//            base = ambiguousCallChar;
        } else if (tt[0] && !tt[1] && !tt[2] && !tt[3]) { //A
            base = 'A';
        } else if (!tt[0] && tt[1] && !tt[2] && !tt[3]) { //C
            base = 'C';
        } else if (!tt[0] && !tt[1] && tt[2] && !tt[3]) { //G
            base = 'G';
        } else if (!tt[0] && !tt[1] && !tt[2] && tt[3]) { //T
            base = 'T';
        } else if (tt[0] && !tt[1] && tt[2] && !tt[3]) { //A or G 
            base = 'R';
        } else if (!tt[0] && tt[1] && !tt[2] && tt[3]) { //C or T 
            base = 'Y';
        } else if (!tt[0] && tt[1] && tt[2] && !tt[3]) { //G or C 
            base = 'S';
        } else if (tt[0] && !tt[1] && !tt[2] && tt[3]) { //A or T 
            base = 'W';
        } else if (!tt[0] && !tt[1] && tt[2] && tt[3]) { //G or T 
            base = 'K';
        } else if (tt[0] && tt[1] && !tt[2] && !tt[3]) { //A or C 
            base = 'M';
        } else if (!tt[0] && tt[1] && tt[2] && tt[3]) { //C or G or T
            base = 'B';
        } else if (tt[0] && !tt[1] && tt[2] && tt[3]) { //A or G or T
            base = 'D';
        } else if (tt[0] && tt[1] && !tt[2] && tt[3]) { //A or C or T
            base = 'H';
        } else if (tt[0] && tt[1] && tt[2] && !tt[3]) { //A or C or G
            base = 'V';
        } else if (tt[0] && tt[1] && tt[2] && tt[3]) {   //any base
            base = 'N';
        }

//        int maxErrLocusInt = (int) Math.floor(locusDepth * maxPercErrLocus / 100);
//        if (totalOther > maxErrLocusInt) {
//            return Character.toLowerCase(base);
//        } else {
            return base;
//        }
//        return '#'; //should never happen
    }
}
