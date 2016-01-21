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

import shared.SequenceOps;

/**
 * One of the two representations of a k-mer, with one of the ends (typically 1
 * base) stored separately from the reminder (core - typically k-1 bases). The
 * object has the core field filled and either the leftCLip or the rightCLip
 * Only used as a wrapper class when generating more compact PairMer
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerGenerator {

//    private final String leftClip;
//    private final String core;
//    private final String rightClip;
    private static final int MAX_1LONG_ENCODE = 32; //2bits per nucl, signed long so should be 31, but can use sign bit if lex ordering not needed, so 32 allowed 
    private static final int MAX_2LONG_ENCODE = 64;
    private static final int MAX_3LONG_ENCODE = 96;
    private static final int MAX_4LONG_ENCODE = 128;
    private static final int MAX_5LONG_ENCODE = 160;

    public static PairMer generatePairMer(String kmerString, boolean frontClip, int overlapLength) {
        
        char leftClip = '#';
        String core;
        char rightClip = '#';

        //SPLIT THE INPUT INTO CORE AND CLIP
        String coreTmp;
        char clip;
        int len = kmerString.length();
        if (frontClip) {
            coreTmp = kmerString.substring(len - overlapLength);
            clip = kmerString.charAt(0);
//            System.err.println("  "+coreTmp+" <-coreTmp");
//            clip = kmerString.substring(0, len - overlapLength);
        } else {
            coreTmp = kmerString.substring(0, overlapLength);
            clip = kmerString.charAt(overlapLength);
//            System.err.println(" "+coreTmp+"  <-coreTmp");
        }

        //ORIENTATE CORE AND CLIP BASED ON LEX ORDER OF CORE AND ITS REV-COMP
        String coreRC = SequenceOps.getReverseComplementString(coreTmp);
        if (coreRC.compareTo(coreTmp) < 0) {
            //REV_COMP = TRUE
            core = coreRC;
            if (frontClip) {
//                leftClip = ' ';
                rightClip = SequenceOps.complement(clip);
//                System.out.println(clip+"-"+rightClip);
//                rightClip = SequenceOps.getReverseComplementString(clip);
            } else {
//                System.out.println(clip+"-"+leftClip);
                leftClip = SequenceOps.complement(clip);
//                rightClip = ' ';
            }
        } else {
            core = coreTmp;
            if (frontClip) {
                leftClip = clip;
//                rightClip = " ";
            } else {
//                leftClip = " ";
                rightClip = clip;
            }
        }
//        if(leftClip == '#')
//            System.err.println(leftClip+""+core+""+rightClip+" <-Generated");
//        else 
//            System.err.println(" "+leftClip+""+core+""+rightClip+" <-Generated");
            
        if (kmerString.length() - 1 <= MAX_1LONG_ENCODE) {
            return new PairMer1LongEncoded(leftClip, core, rightClip);
        } else if (kmerString.length() - 1 <= MAX_2LONG_ENCODE) {
            return new PairMer2LongEncoded(leftClip, core, rightClip);
        } else if (kmerString.length() - 1 <= MAX_3LONG_ENCODE) {
            return new PairMer3LongEncoded(leftClip, core, rightClip);
        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
            return new PairMer4LongEncoded(leftClip, core, rightClip);
        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
            return new PairMer5LongEncoded(leftClip, core, rightClip);
        } else {
            return new PairMerIntArrEncoded(leftClip, core, rightClip);
        }
    }
    
    /**
     * Given a core string, generates a PairMer (without clips)
     * for querying populated PairMerMaps
     *
     * @param core, which will be converted to its canonical form
     * @param k, used to encode pairMer using the correct data structure
     * @return PairMer to be used for interrogating a Map
     */
    public static PairMer getPairMer(String core, int k) {
        if (k - 1 <= MAX_1LONG_ENCODE) {
            return new PairMer1LongEncoded(core);
        } else if (k - 1 <= MAX_2LONG_ENCODE) {
            return new PairMer2LongEncoded(core);
        } else if (k - 1 <= MAX_3LONG_ENCODE) {
            return new PairMer3LongEncoded(core);
        } else if (k - 1 <= MAX_4LONG_ENCODE) {
            return new PairMer4LongEncoded(core);
        } else if (k - 1 <= MAX_5LONG_ENCODE) {
            return new PairMer5LongEncoded(core);
        } else {
            return new PairMerIntArrEncoded(core);
        }
    }

//    public String getLeftClip() {
//        return leftClip;
//    }
//
//    public char getLeftClipChar() {
//        return leftClip.charAt(0);
//    }
//
//    public String getCore() {
//        return core;
//    }
//
//    public String getRightClip() {
//        return rightClip;
//    }
//
//    public char getRightClipChar() {
//        return rightClip.charAt(0);
//    }
//    public PairMer generatepairMer(String kmerString) {
//        if (kmerString.length() - 1 <= MAX_1LONG_ENCODE) {
//            return new PairMer1LongEncoded();
////        } else if (kmerString.length() - 1 <= MAX_2LONG_ENCODE) {
////            return new PairMer2LongEncoded(this);
////        } else if (kmerString.length() - 1 <= MAX_3LONG_ENCODE) {
////            return new PairMer3LongEncoded(this);
////        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
////            return new PairMer4LongEncoded(this);
////        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
////            return new PairMer5LongEncoded(this);
////        } else {
////            return new PairMerIntArrEncoded(this);
//        }
//    }
}
