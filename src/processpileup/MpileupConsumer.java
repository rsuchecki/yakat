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
import java.util.regex.Pattern;
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
    private final double minMinorMajorRatio;
    private final boolean allWithinThresholds;
//    private final boolean reportHets;
    private final String TOOL_NAME;
    private final PrintStream bufferedOut;

    private final char zeroReadsChar;
    private final char ambiguousCallChar;

    private final Integer minSnpsToRef;

    private final int minCallsHet;
    private final int maxCallsHet;
    private final int minCallsUncertain;
    private final int maxCallsUnceratin;

    private int minMissingSamples;
    private int maxUnCalledSamples;

    public MpileupConsumer(BlockingQueue<ArrayList<String>> inputQueue, OptSet optSet, String TOOL_NAME, PrintStream bufferedOut) {
        minCoveragePerLocus = (int) optSet.getOpt("c").getValueOrDefault();
        maxCoveragePerLocus = (int) optSet.getOpt("C").getValueOrDefault();
        minSamples = (int) optSet.getOpt("s").getValueOrDefault();
        minCoveragePerAllele = (int) optSet.getOpt("a").getValueOrDefault();
        maxPercErrAllele = (double) optSet.getOpt("A").getValueOrDefault();
        maxPercErrLocus = (double) optSet.getOpt("L").getValueOrDefault();
//        OUT_BUFFER_SIZE = (int) optSet.getOpt("u").getValueOrDefault();
        minMinorMajorRatio = (double) optSet.getOpt("min-minor-major-ratio").getValueOrDefault();

        zeroReadsChar = (Character) optSet.getOpt("zero-reads-char").getValueOrDefault();
        ambiguousCallChar = (Character) optSet.getOpt("ambiguous-call-char").getValueOrDefault();

        minSnpsToRef = (Integer) optSet.getOpt("min-snps-to-reference").getValueOrDefault();

        minCallsHet = (int) optSet.getOpt("min-calls-het").getValueOrDefault();
        maxCallsHet = (int) optSet.getOpt("max-calls-het").getValueOrDefault();
        minCallsUncertain = (int) optSet.getOpt("min-calls-uncertain").getValueOrDefault();
        maxCallsUnceratin = (int) optSet.getOpt("max-calls-uncertain").getValueOrDefault();

        allWithinThresholds = optSet.getOptS('W').isUsed();
//        reportHets = optSet.getOptS('H').isUsed();
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
            Pattern spliPattern = Pattern.compile("\t");
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    String[] toks = spliPattern.split(line);
                    if (toks.length == 4 && toks[3].equals("0")) {
                        continue;
                    }
//                    String[] toks = line.split("\t");
                    //toks 0,1,2 are ref, position, refbase
                    char refBase = toks[2].charAt(0);
                    char refBaseUpper = Character.toUpperCase(refBase);
                    StringBuilder coveragesSB = new StringBuilder("COUNTS");
                    StringBuilder callsSB = new StringBuilder("CALLS");
                    StringBuilder common = new StringBuilder("\t");
                    common.append(toks[0]).append("\t").append(toks[1]).append("\t").append(toks[2]);
                    callsSB.append(common);
                    coveragesSB.append(common);
                    boolean betweenSamplesSnps = false;
//                    String lastCalledBase = "N";
                    char lastCalledBase = '#';
//                    int samplesWithBaseCalled = 0;
//                    int[] basesAllSamples = new int[7];

                    int coverageAllSamples = 0;
                    int samplesWithinCoverage = 0;
                    int uncertain = 0;
                    int samplesZeroCoverage = 0;
                    int snpsToRef = 0;
                    int hetCalls = 0;
                    for (int i = 3; i < toks.length; i += 3) {
                        try {
                            coveragesSB.append("\t");
                            callsSB.append("\t");
                            int coverage = Integer.parseInt(toks[i]);
                            if (coverage == 0) {
                                coveragesSB.append("0,0,0,0,0,0");
                                callsSB.append(zeroReadsChar);
                                continue;
                            }

                            coverageAllSamples += coverage;
                            if (coverage >= minCoveragePerLocus && coverage <= maxCoveragePerLocus) {
                                samplesWithinCoverage++;
                            }
                            int[] bases = getBaseCounts(toks[i + 1], refBase);
//                            for (int j = 0; j < bases.length; j++) {
//                                basesAllSamples[j] += bases[j];
//                            }

                            char calledBase = getIUPAC(bases, coverage);
                            char calledBaseUpper = Character.toUpperCase(calledBase);

                            if (calledBase == ambiguousCallChar) {
                                uncertain++;
                            } else if (calledBase == zeroReadsChar) {
                                samplesZeroCoverage++;
                            } else {
                                if (!isACGT(calledBase)) {
                                    hetCalls++;
                                }
                                if (calledBaseUpper != refBaseUpper) {
                                    snpsToRef++;
                                }
                                if (lastCalledBase != '#' && Character.toUpperCase(calledBase) != Character.toUpperCase(lastCalledBase)) {
                                    betweenSamplesSnps = true;
                                }
                                lastCalledBase = calledBase;
//                                samplesWithBaseCalled++;
                            }
                            callsSB.append(calledBase);

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
//                    if (samplesWithinCoverage >= minSamples && (calledDifferentBases || (allWithinThresholds && basesCalled>0 && unknown==0 && uncovered>1))) {
                    if (samplesWithinCoverage >= minSamples) {
                        if (betweenSamplesSnps || (minSnpsToRef != null && snpsToRef >= minSnpsToRef)
                                || (allWithinThresholds
                                && hetCalls >= minCallsHet && hetCalls <= maxCallsHet
                                && uncertain >= minCallsUncertain && uncertain <= maxCallsUnceratin)) {
//                            && samplesWithBaseCalled > 0 && uncertain == 0 && samplesZeroCoverage > 1)) {
//                            coveragesSB.append("\t");
//                            for (int j = 1; j < basesAllSamples.length; j++) {
//                                coveragesSB.append(basesAllSamples[j]);
//                                if (j < basesAllSamples.length - 1) {
//                                    coveragesSB.append(",");
//                                }
//                            }
//                            callsSB.append("\t").append(getIUPAC(basesAllSamples, coverageAllSamples));

                            bufferedOut.println(coveragesSB);
                            bufferedOut.println(callsSB);
//                        System.out.println(coveragesSB);
//                        System.out.println(callsSB);
//                        bufferList.add(coveragesSB.toString());
//                        bufferList.add(callsSB.toString());
//                        if (bufferList.size() == OUT_BUFFER_SIZE) {
//                            outputQueue.put(bufferList);
//                            bufferList = new ArrayList<>();
                        }
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

    private boolean isACGT(char c) {
        switch (c) {
            case 'A':
                return true;
            case 'C':
                return true;
            case 'G':
                return true;
            case 'T':
                return true;
            case 'a':
                return true;
            case 'c':
                return true;
            case 'g':
                return true;
            case 't':
                return true;
            default:
                return false;
        }
    }

    /**
     * Given an array of numbers representing observed base counts
     * [?,#A,#C,#G,#T] (index zero ignored for now) call IUPAC representation of
     * the observed bases. Bases are not considered for IUPAC reporting if there
     * are less than minCoverageThreshold or
     *
     * @param encoded [?,#A,#C,#G,#T] (index zero ignored for now)
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
    private char getIUPAC(int[] encoded, int locusDepth) {
        if (locusDepth == 0) {
            return zeroReadsChar;
        } else if (locusDepth < minCoveragePerLocus || locusDepth > maxCoveragePerLocus) {
            return ambiguousCallChar;
        }
        boolean tt[] = new boolean[7]; //truth table i=1,2,3,4 values A,C,G,T...
        //Convert max perc error into max int coverage
        int maxErrAlleleInt = (int) Math.floor(locusDepth * maxPercErrAllele / 100);
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
        char base = ambiguousCallChar;
        if (totalOther > maxErrAlleleInt) {
            return ambiguousCallChar;
        }
        if (!tt[1] && !tt[2] && !tt[3] && !tt[4]) {
            if (tt[5] && !tt[6]) { //insertion in reference
                base = 'E';
            } else if (!tt[5] && tt[6]) { //deletion in reference
                base = 'I';
            } else if (tt[5] && tt[6]) { //deletion and insertion??
                base = ambiguousCallChar;
            } else if (totalOther == 0) {
                base = zeroReadsChar;
            } else {
                base = ambiguousCallChar;
            }
        } else if ((tt[1] || tt[2] || tt[3] || tt[4]) && (tt[5] || tt[6])) { //base and an indel???
            base = ambiguousCallChar;
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

        int maxErrLocusInt = (int) Math.floor(locusDepth * maxPercErrLocus / 100);
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

//        boolean processIndel = false;
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '.':
                case ',':
                    c = refBase;
            }

            //Not processing an indel
            if (lenInsRef == null && lenDelRef == null) {
//            if (!processIndel) {
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
                    case '+': //insertion in the read indicator, pattern: `\+[0-9]+[ACGTNacgtn]+'            
                        lenInsRef = new StringBuilder();
//                        processIndel = true;
                        break;
                    case '-': //deletion  in the read indicator, pattern: `\-[0-9]+[ACGTNacgtn]+'            
                        lenDelRef = new StringBuilder();
//                        processIndel = true;
                        break;
                    default:
//                    bases[0]++;
                }
            } else {//processing an indel                        
                if (lenInsRef != null) {
                    if (Character.isDigit(c)) {
                        lenInsRef.append(c);
                    } else {
//                        switch (c) {
//                            case 'A':
//                            case 'a':
//                                break;
//                            case 'C':
//                            case 'c':
//                                break;
//                            case 'G':
//                            case 'g':
//                                break;
//                            case 'T':
//                            case 't':
//                                break;
//                            case 'N':
//                            case 'n':
//                                break;
//                        }
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
