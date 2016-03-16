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
package processpileup;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class MpileupConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final int minCoverageThreshold;
    private final int maxCoverageThreshold;
    private final int minSamples;
    private final int maxPercAlternative;
    private final String TOOL_NAME;

    public MpileupConsumer(BlockingQueue<ArrayList<String>> inputQueue, int minCoverageThreshold,
            int maxCoverageThreshold, int minSamples, int maxPercAlternative, String TOOL_NAME) {
        this.inputQueue = inputQueue;
        this.minCoverageThreshold = minCoverageThreshold;
        this.maxCoverageThreshold = maxCoverageThreshold;
        this.minSamples = minSamples;
        this.maxPercAlternative = maxPercAlternative;
        this.TOOL_NAME = TOOL_NAME;
    }

    @Override
    public void run() {
        try {
            ArrayList<String> list;
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    String[] toks = line.split("\t");
                    int r = 0;
                    //tos 0,1,2 are ref, position, refbase
                    char refBase = toks[2].charAt(0);
                    StringBuilder coveragesSB = new StringBuilder("COUNTS");
                    StringBuilder callsSB = new StringBuilder("CALLS");
                    StringBuilder common = new StringBuilder("\t");
                    common.append(toks[0]).append("\t").append(toks[1]).append("\t").append(toks[2]);
                    callsSB.append(common);
                    coveragesSB.append(common);
                    boolean calledDifferentBases = false;
//                    String lastCalledBase = "N";
                    char lastCalledBase = '?';

                    int samplesWithinCoverage = 0;
                    for (int i = 3; i < toks.length; i += 3) {
                        try {
                            int coverage = Integer.parseInt(toks[i]);
                            if (coverage >= minCoverageThreshold && coverage <= maxCoverageThreshold) {
                                samplesWithinCoverage++;
                            }
                            int[] bases = getBaseCounts(toks[i + 1], refBase);
                            coveragesSB.append("\t");
//                            String calledBase = callBase(bases, maxPercAlternative, minCoverageThreshold);

                            char calledBase = getIUPAC(bases, coverage, maxPercAlternative, minCoverageThreshold);
//                            if (calledBase.matches("A|T|C|G")) {
                            if (calledBase != '?' && calledBase != ' ') {
//                                if (!calledBase.equals(lastCalledBase) && lastCalledBase.matches("A|T|C|G")) {
                                if (calledBase != lastCalledBase && lastCalledBase != '?' && lastCalledBase != ' ') {
                                    calledDifferentBases = true;
                                }
                                lastCalledBase = calledBase;
                            }
                            callsSB.append("\t").append(calledBase);
                            for (int j = 1; j < bases.length; j++) {
                                coveragesSB.append(bases[j]);
                                if (j < bases.length - 1) {
                                    coveragesSB.append(",");
                                }
                            }
                        } catch (ArrayIndexOutOfBoundsException e) {
                            Reporter.report("[FATAL]", "Array index out of bounds - likely cause: mismatch between samples given and pileup file content", TOOL_NAME);
                            System.exit(1);
                        }
                    }
                    if (samplesWithinCoverage >= minSamples && calledDifferentBases) {
                        System.out.println(coveragesSB);
                        System.out.println(callsSB);
                    }
                }
            }
            inputQueue.put(new ArrayList<String>()); //inform other threads
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

    }

    private char callBase(int[] bases, int maxPercAlternativeAllowed, int minCoverageThreshold) {
        int maxCov = 0;
        int maxPos = 0;
        int totalDepth = 0;
        //Identify "major allele"
        for (int i = 1; i < bases.length; i++) {
            if (bases[i] > maxCov) {
                maxCov = bases[i];
                maxPos = i;
            }
            totalDepth += bases[i];
        }

        double percentOfOtherBases = 1 - ((double) maxCov / totalDepth);

        if (maxCov >= minCoverageThreshold) {
            if (percentOfOtherBases > maxPercAlternativeAllowed / 100) {
                return 'N'; //
            }
            return getBase(maxPos);
        }
        if (maxCov == 0) {
            return ' ';
        }
        //If uncertain
        return '?';
    }

//    private char getIUPAC(int [] callCounts, int ) {
//        
//    }
    private char getBase(int encoded) {
        switch (encoded) {
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 3:
                return 'G';
            case 4:
                return 'T';
            default:
                return 'N';
        }
    }

    /**
     * Given an array of numbers representing observed base counts  [?,#A,#C,#G,#T] (index zero ignored for now) 
     * call IUPAC representation of the observed bases. Bases are not considered for IUPAC reporting if there are less 
     * than minCoverageThreshold or 
     * @param encoded [?,#A,#C,#G,#T] (index zero ignored for now)
     * @param totalDepth 
     * @param maxErrorPercent % of totalDepth threshold above which a base 
     * (or resulting IUPAC) is reported rather than being ignored as erroneous
     * @param minCoverageThreshold int min depth for a base to be considered
     * 
     * TODO EXTEND TO ACCOMMODATE IDELS
     * 
     * @return 
     */
    private char getIUPAC(int[] encoded, int totalDepth, int maxErrorPercent, int minCoverageThreshold) {
        boolean tt[] = new boolean[5]; //truth table i=1,2,3,4 values A,C,G,T...
        //Convert max perc error into max int coverage
        int maxErrInt = (int) Math.floor(totalDepth * (float) maxErrorPercent / 100);        
        //Max of maxErr min coverage threshold 
//        int minCoverage = Math.max((int) Math.ceil(totalDepth * (float) maxErrorPercent / 100), minCoverageThreshold);
        int totalOther = 0;
        for (int i = 1; i < encoded.length; i++) {
            if (encoded[i] >= minCoverageThreshold && encoded[i] > maxErrInt) {
                tt[i] = true;                
            } else {
                totalOther += encoded[i];
            }
        }
        if (totalOther > maxErrInt) {
            return '?';
        }
        if (!tt[1] && !tt[2] && !tt[3] && !tt[4]) {            
            if(totalOther == 0)
                return ' ';
            else 
                return  '?';
        } else if (tt[1] && !tt[2] && !tt[3] && !tt[4]) { //A
            return 'A';
        } else if (!tt[1] && tt[2] && !tt[3] && !tt[4]) { //C
            return 'C';
        } else if (!tt[1] && !tt[2] && tt[3] && !tt[4]) { //G
            return 'G';
        } else if (!tt[1] && !tt[2] && !tt[3] && tt[4]) { //T
            return 'T';
        } else if (tt[1] && !tt[2] && tt[3] && !tt[4]) { //A or G 
            return 'R';
        } else if (!tt[1] && tt[2] && !tt[3] && tt[4]) { //C or T 
            return 'Y';
        } else if (!tt[1] && tt[2] && tt[3] && !tt[4]) { //G or C 
            return 'S';
        } else if (tt[1] && !tt[2] && !tt[3] && tt[4]) { //A or T 
            return 'W';
        } else if (!tt[1] && !tt[2] && tt[3] && tt[4]) { //G or T 
            return 'K';
        } else if (tt[1] && tt[2] && !tt[3] && !tt[4]) { //A or C 
            return 'M';
        } else if (!tt[1] && tt[2] && tt[3] && tt[4]) { //C or G or T
            return 'B';
        } else if (tt[1] && !tt[2] && tt[3] && tt[4]) { //A or G or T
            return 'D';
        } else if (tt[1] && tt[2] && !tt[3] && tt[4]) { //A or C or T
            return 'H';
        } else if (tt[1] && tt[2] && tt[3] && !tt[4]) { //A or C or G
            return 'V';
        } else if (tt[1] && tt[2] && tt[3] && tt[4]) {   //any base
            return 'N';
        }
        return '#'; //should never happen
    }

    private int[] getBaseCounts(String s, char refBase) {
//        s = s.replaceAll("\\.|,", String.valueOf(refBase)).replaceAll("\\*", "I").replaceAll("\\^.|\\$", "");
        int[] bases = new int[5]; // ?,A,C,G,T,insRef,delRef 
        StringBuilder lenInsRef = null;
        StringBuilder lenDelRef = null;

        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '.':
                case ',':
                    c = refBase;
            }

            //Not processing an indel
            if (lenInsRef == null && lenDelRef == null) {
                switch (c) {
                    case 'A':
                    case 'a':
                        bases[1]++;
                        break;
                    case 'C':
                    case 'c':
                        bases[2]++;
                        break;
                    case 'G':
                    case 'g':
                        bases[3]++;
                        break;
                    case 'T':
                    case 't':
                        bases[4]++;
                        break;
                    case '^': //aligned fragment start
                    case '$': //aligned fragment end
                        break;
                    case '+': //insertion in the reference indicator, pattern: `\+[0-9]+[ACGTNacgtn]+'            
                        lenInsRef = new StringBuilder();
                        break;
                    case '-': //deletion  in the reference indicator, pattern: `\-[0-9]+[ACGTNacgtn]+'            
                        lenDelRef = new StringBuilder();
                        break;
                    default:
//                    bases[0]++;
                }
            } else {//processing an indel                        
                if (lenInsRef != null) {
                    if (Character.isDigit(c)) {
                        lenInsRef.append(c);
                    } else {
                        int indelLen = Integer.parseInt(lenInsRef.toString());
                        lenInsRef = null;
                        i += indelLen - 1;
                        bases[5]++;
                    }
                }
                if (lenDelRef != null) {
                    if (Character.isDigit(c)) {
                        lenDelRef.append(c);
                    } else {
                        int indelLen = Integer.parseInt(lenDelRef.toString());
                        lenDelRef = null;
                        i += indelLen - 1;
                        bases[6]++;
                    }
                }
            }
        }
        return bases;
    }
}
