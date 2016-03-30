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

import argparser.OptSet;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class MpileupConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
//    private final BlockingQueue<ArrayList<String>> outputQueue;
//    private final int OUT_BUFFER_SIZE;
    private final int minCoveragePerLocus;
    private final int maxCoveragePerLocus;
    private final int minSamples;
    private final int minCoveragePerAllele;
    private final double maxPercErrAllele;
    private final double maxPercErrLocus;
    private final boolean allWithinThresholds;
    private final String TOOL_NAME;
    private final PrintStream bufferedOut;

    public MpileupConsumer(BlockingQueue<ArrayList<String>> inputQueue, OptSet optSet, String TOOL_NAME, PrintStream bufferedOut) {
        minCoveragePerLocus = (int) optSet.getOpt("c").getValueOrDefault();
        maxCoveragePerLocus = (int) optSet.getOpt("C").getValueOrDefault();
        minSamples = (int) optSet.getOpt("s").getValueOrDefault();
        minCoveragePerAllele = (int) optSet.getOpt("a").getValueOrDefault();
        maxPercErrAllele = (double) optSet.getOpt("A").getValueOrDefault();
        maxPercErrLocus = (double) optSet.getOpt("L").getValueOrDefault();
//        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
        allWithinThresholds = optSet.getOptS('W').isUsed();
        this.inputQueue = inputQueue;
//        this.outputQueue = outputQueue;
        this.TOOL_NAME = TOOL_NAME;
        this.bufferedOut = bufferedOut;
    }

    
//    public MpileupConsumer(BlockingQueue<ArrayList<String>> inputQueue, int minCoveragePerLocus,
//        int maxCoverageThreshold, int minSamples, double maxPercAlternative, String TOOL_NAME) {
//        this.inputQueue = inputQueue;
//        
////                int minCoveragePerLocus = (int) optSet.getOpt("c").getValueOrDefault();
////        int maxCoveragePerLocus = (int) optSet.getOpt("C").getValueOrDefault();
////        int minSamples = (int) optSet.getOpt("s").getValueOrDefault();
////        int minCoveragePerAllele = (int) optSet.getOpt("a").getValueOrDefault();
////        
////        double maxPercErrAllele = (double) optSet.getOpt("A").getValueOrDefault();
////        double maxPercErrLocus = (double) optSet.getOpt("L").getValueOrDefault();
//        
//        
//        this.minCoverageThreshold = minCoverageThreshold;
//        this.maxCoverageThreshold = maxCoverageThreshold;
//        this.minSamples = minSamples;
//        this.maxPercAlternative = maxPercAlternative;
//        this.TOOL_NAME = TOOL_NAME;
//    }
    @Override
    public void run() {
        try {
//            ArrayList<String> bufferList = new ArrayList<>(OUT_BUFFER_SIZE);
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
                    char lastCalledBase = '#';
                    int basesCalled = 0;
                    int[] basesAllSamples = new int[7];

                    int coverageAllSamples = 0;
                    int samplesWithinCoverage = 0;
                    int unknown = 0;
                    int uncovered = 0;
                    for (int i = 3; i < toks.length; i += 3) {
                        try {
                            int coverage = Integer.parseInt(toks[i]);
                            coverageAllSamples += coverage;
                            if (coverage >= minCoveragePerLocus && coverage <= maxCoveragePerLocus) {
                                samplesWithinCoverage++;
                            }
                            int[] bases = getBaseCounts(toks[i + 1], refBase);
                            for (int j = 0; j < bases.length; j++) {
                                basesAllSamples[j] += bases[j];
                            }
                            coveragesSB.append("\t");
//                            String calledBase = callBase(bases, maxPercAlternative, minCoverageThreshold);

                            char calledBase = getIUPAC(bases, coverage);
//                            if (calledBase.matches("A|T|C|G")) {
                            if(calledBase == '?') {
                                unknown++;
                            } else if (calledBase == '.') {
                                uncovered++;
                            } else {
                                if (Character.toUpperCase(calledBase) != Character.toUpperCase(lastCalledBase) && lastCalledBase != '#') {
                                    calledDifferentBases = true;
                                }
                                lastCalledBase = calledBase;
                                basesCalled++;
                            }
                            callsSB.append("\t");
                            callsSB.append(calledBase);
//                            if (calledBase != ' ') { //' ' != ""   which may have an impact on downstream analysis
//                            }
                            for (int j = 1; j < bases.length; j++) {
                                coveragesSB.append(bases[j]);
                                if (j < bases.length - 1) {
                                    coveragesSB.append(",");
                                }
                            }
//                            coveragesSB.append(System.lineSeparator());
//                            callsSB.append(System.lineSeparator());
                        } catch (ArrayIndexOutOfBoundsException e) {
                            Reporter.report("[FATAL]", "Array index out of bounds - likely cause: mismatch between samples given and pileup file content", TOOL_NAME);
                            System.err.println(line);

                            System.exit(1);
                        }
                    }
                    //TODO uncovered should be set to > 0, > 1 is specific to having a very poor sample here
                    if (samplesWithinCoverage >= minSamples && (calledDifferentBases || (allWithinThresholds && basesCalled>0 && unknown==0 && uncovered>1))) {

                        coveragesSB.append("\t");
                        for (int j = 1; j < basesAllSamples.length; j++) {
                            coveragesSB.append(basesAllSamples[j]);
                            if (j < basesAllSamples.length - 1) {
                                coveragesSB.append(",");
                            }
                        }
                        callsSB.append("\t").append(getIUPAC(basesAllSamples, coverageAllSamples));

                        bufferedOut.println(coveragesSB);
                        bufferedOut.println(callsSB);
//                        System.out.println(coveragesSB);
//                        System.out.println(callsSB);
//                        bufferList.add(coveragesSB.toString());
//                        bufferList.add(callsSB.toString());
//                        if (bufferList.size() == OUT_BUFFER_SIZE) {
//                            outputQueue.put(bufferList);
//                            bufferList = new ArrayList<>();
//                        }
                    }
                }
                bufferedOut.flush();
            }
//            outputQueue.put(bufferList);
//            outputQueue.put(new ArrayList<String>()); //inform other threads
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
     * Given an array of numbers representing observed base counts [?,#A,#C,#G,#T] (index zero ignored for now) call
     * IUPAC representation of the observed bases. Bases are not considered for IUPAC reporting if there are less than
     * minCoverageThreshold or
     *
     * @param encoded [?,#A,#C,#G,#T] (index zero ignored for now)
     * @param totalDepth
     * @param maxErrorPercent % of totalDepth threshold above which a base (or resulting IUPAC) is reported rather than
     * being ignored as erroneous
     * @paramminCoveragePerAlleleminCoverageThreshold int min depth for a base to be considered
     *
     * TODO EXTEND TO ACCOMMODATE IDELS
     *
     * @return
     */
    private char getIUPAC(int[] encoded, int totalDepth) {
        boolean tt[] = new boolean[7]; //truth table i=1,2,3,4 values A,C,G,T...
        //Convert max perc error into max int coverage
        int maxErrAlleleInt = (int) Math.floor(totalDepth * maxPercErrAllele / 100);
//        System.err.println("MaxErrInt="+maxErrInt);
        //Max of maxErr min coverage threshold 
//        int minCoverage = Math.max((int) Math.ceil(totalDepth * (float) maxErrorPercent / 100), minCoverageThreshold);
        int totalOther = 0;
        for (int i = 1; i < encoded.length; i++) {
//            if (encoded[i] > 0) {
//                tt[i] = true;
//            }
//            if (encoded[i] < minCoverageThreshold || encoded[i] <= maxErrInt) {
//                totalOther += encoded[i];                
//            }

            if (encoded[i] >= minCoveragePerAllele && encoded[i] > maxErrAlleleInt) {
                tt[i] = true;
            } else {
                totalOther += encoded[i];
            }
        }
        char base = '?';
        if (totalOther > maxErrAlleleInt) {
            return '?';
        }
        if (!tt[1] && !tt[2] && !tt[3] && !tt[4]) {
            if (tt[5] && !tt[6]) { //insertion in reference
                base = 'E';
            } else if (!tt[5] && tt[6]) { //deletion in reference
                base = 'I';
            } else if (tt[5] && tt[6]) { //deletion and insertion??
                base = '?';
            } else if (totalOther == 0) {
                base = '.';
            } else {
                base = '?';
            }
        } else if ((tt[1] || tt[2] || tt[3] || tt[4]) && (tt[5] || tt[6])) { //base and an indel???
            base = '?';
        } else if (tt[1] && !tt[2] && !tt[3] && !tt[4]) { //A
            base = 'A';
        } else if (!tt[1] && tt[2] && !tt[3] && !tt[4]) { //C
            base = 'C';
        } else if (!tt[1] && !tt[2] && tt[3] && !tt[4]) { //G
            base = 'G';
        } else if (!tt[1] && !tt[2] && !tt[3] && tt[4]) { //T
            base = 'T';
        } else if (tt[1] && !tt[2] && tt[3] && !tt[4]) { //A or G 
            base = 'R';
        } else if (!tt[1] && tt[2] && !tt[3] && tt[4]) { //C or T 
            base = 'Y';
        } else if (!tt[1] && tt[2] && tt[3] && !tt[4]) { //G or C 
            base = 'S';
        } else if (tt[1] && !tt[2] && !tt[3] && tt[4]) { //A or T 
            base = 'W';
        } else if (!tt[1] && !tt[2] && tt[3] && tt[4]) { //G or T 
            base = 'K';
        } else if (tt[1] && tt[2] && !tt[3] && !tt[4]) { //A or C 
            base = 'M';
        } else if (!tt[1] && tt[2] && tt[3] && tt[4]) { //C or G or T
            base = 'B';
        } else if (tt[1] && !tt[2] && tt[3] && tt[4]) { //A or G or T
            base = 'D';
        } else if (tt[1] && tt[2] && !tt[3] && tt[4]) { //A or C or T
            base = 'H';
        } else if (tt[1] && tt[2] && tt[3] && !tt[4]) { //A or C or G
            base = 'V';
        } else if (tt[1] && tt[2] && tt[3] && tt[4]) {   //any base
            base = 'N';
        }

        int maxErrLocusInt = (int) Math.floor(totalDepth * maxPercErrLocus / 100);
        if (totalOther > maxErrLocusInt) {
            return Character.toLowerCase(base);
        } else {
            return base;
        }
//        return '#'; //should never happen
    }

    private int[] getBaseCounts(String s, char refBase) {
//        s = s.replaceAll("\\.|,", String.valueOf(refBase)).replaceAll("\\*", "I").replaceAll("\\^.|\\$", "");
        int[] bases = new int[7]; // ?,A,C,G,T,insRef,delRef 
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
//                        System.out.println("Ins in ref "+s);
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
//                        System.out.println("del in ref"+s);
                    }
                }
            }
        }
        return bases;
    }
}
