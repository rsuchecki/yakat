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
    private static final int MAX_1LONG_ENCODE = 32; //using all bits, previously: 2bits per nucl, signed long so should be 31, but can use sign bit if lex ordering not needed
    private static final int MAX_2LONG_ENCODE = 64;
    private static final int MAX_3LONG_ENCODE = 96;
//    private static final int MAX_4LONG_ENCODE = 128 - 1;
//    private static final int MAX_5LONG_ENCODE = 160 - 1;

//    /**
//     * GeneRatePairMer object with k-1 core separated from 1 base clip 
//     * Clip end is selected based on frontClip param
//     * 
//     * @param kmerString
//     * @param frontClip 
//     * @param overlapLength - i.e. k-1
//     * @return 
//     */
//    public static PairMer generatePairMer(String kmerString, boolean frontClip, int overlapLength) {        
//        char leftClip = '#'; 
//        String core;
//        char rightClip = '#';
//        
//        
//
//        //SPLIT THE INPUT INTO CORE AND CLIP
//        String coreTmp;
//        char clip;
//        int len = kmerString.length();
//        if (frontClip) {
//            coreTmp = kmerString.substring(len - overlapLength);
//            clip = kmerString.charAt(0);
//        } else {
//            coreTmp = kmerString.substring(0, overlapLength);
//            clip = kmerString.charAt(overlapLength);
//        }
//
//        //ORIENTATE CORE AND CLIP BASED ON LEX ORDER OF CORE AND ITS REV-COMP
//        String coreRC = SequenceOps.getReverseComplementString(coreTmp);
//        if (coreRC.compareTo(coreTmp) < 0) {
//            //REV_COMP = TRUE
//            core = coreRC;
//            if (frontClip) {
//                rightClip = SequenceOps.complement(clip);
//            } else {
//                leftClip = SequenceOps.complement(clip);
//            }
//        } else {
//            core = coreTmp;
//            if (frontClip) {
//                leftClip = clip;
//            } else {
//                rightClip = clip;
//            }
//        }
//            
////        if (kmerString.length() - 1 <= MAX_1LONG_ENCODE) {
////            return new PairMer1LongEncoded(leftClip, core, rightClip);
////        } else if (kmerString.length() - 1 <= MAX_2LONG_ENCODE) {
////            return new PairMer2LongEncoded(leftClip, core, rightClip);
////        } else if (kmerString.length() - 1 <= MAX_3LONG_ENCODE) {
////            return new PairMer3LongEncoded(leftClip, core, rightClip);
////        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
////            return new PairMer4LongEncoded(leftClip, core, rightClip);
////        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
////            return new PairMer5LongEncoded(leftClip, core, rightClip);
////        } else {
//            return new PairMerIntArrEncoded(leftClip, core, rightClip);
////        }
//    }
    /**
     * GeneratePairMer object with k-1 core separated from 1 base clip
     *
     * @param sequence
     * @param kmerFrom
     * @param kmerTo inclusive
     * @param frontClip
     * @param freq
     * @return
     * @throws kmerextender.NonACGTException
     */
    public static PairMer generatePairMer(CharSequence sequence, int kmerFrom, int kmerTo, boolean frontClip, int freq) throws NonACGTException {

        if (kmerTo -kmerFrom <= MAX_1LONG_ENCODE) {
            return new PairMer1LongEncoded(sequence, kmerFrom, kmerTo, frontClip, freq);
        } else if (kmerTo -kmerFrom <= MAX_2LONG_ENCODE) {
            return new PairMer2LongEncoded(sequence, kmerFrom, kmerTo, frontClip, freq);
        } else if (kmerTo -kmerFrom  <= MAX_3LONG_ENCODE) { // && kmerTo -kmerFrom > MAX_2LONG_ENCODE) {
            return new PairMer3LongEncoded(sequence, kmerFrom, kmerTo, frontClip, freq);
//        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
//            return new PairMer4LongEncoded(leftClip, core, rightClip);
//        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
//            return new PairMer5LongEncoded(leftClip, core, rightClip);
        } else {
            return new PairMerIntArrEncoded(sequence, kmerFrom, kmerTo, frontClip, freq);
        }
    }

    /**
     * Given a core string, generates a PairMer (without clips) for querying
     * populated PairMerMaps
     *
     * @param core, which will be converted to its canonical form
     * @param k, used to encode pairMer using the correct data structure
     * @return PairMer to be used for interrogating a Map
     * @throws kmerextender.NonACGTException
     */
    public static PairMer getPairMer(CharSequence core, int k) throws NonACGTException {
        if (k - 1 <= MAX_1LONG_ENCODE) {
            return new PairMer1LongEncoded(core);
        } else if (k - 1 <= MAX_2LONG_ENCODE) {
            return new PairMer2LongEncoded(core);
        } else if (k - 1 <= MAX_3LONG_ENCODE ) { // && k-1 > MAX_2LONG_ENCODE) {
            return new PairMer3LongEncoded(core);
//        } else if (k - 1 <= MAX_4LONG_ENCODE) {
//            return new PairMer4LongEncoded(core);
//        } else if (k - 1 <= MAX_5LONG_ENCODE) {
//            return new PairMer5LongEncoded(core);
        } else {
            return new PairMerIntArrEncoded(core);
        }
    }

//    public static PairMer generatePairMer(CharSequence kmerContainingString, int from , int to, boolean frontClip, int overlapLength) {        
//        char leftClip = '#';
//        String core;
//        char rightClip = '#';
//
//        //SPLIT THE INPUT INTO CORE AND CLIP
//        String coreTmp;
//        char clip;
//        int len = to - from +1;
//        if (frontClip) {
//            coreTmp = kmerString.substring(len - overlapLength);
//            clip = kmerString.charAt(0);
//        } else {
//            coreTmp = kmerString.substring(0, overlapLength);
//            clip = kmerString.charAt(overlapLength);
//        }
//
//        //ORIENTATE CORE AND CLIP BASED ON LEX ORDER OF CORE AND ITS REV-COMP
//        String coreRC = SequenceOps.getReverseComplementString(coreTmp);
//        if (coreRC.compareTo(coreTmp) < 0) {
//            //REV_COMP = TRUE
//            core = coreRC;
//            if (frontClip) {
//                rightClip = SequenceOps.complement(clip);
//            } else {
//                leftClip = SequenceOps.complement(clip);
//            }
//        } else {
//            core = coreTmp;
//            if (frontClip) {
//                leftClip = clip;
//            } else {
//                rightClip = clip;
//            }
//        }
//            
//        if (kmerString.length() - 1 <= MAX_1LONG_ENCODE) {
//            return new PairMer1LongEncoded(leftClip, core, rightClip);
//        } else if (kmerString.length() - 1 <= MAX_2LONG_ENCODE) {
//            return new PairMer2LongEncoded(leftClip, core, rightClip);
//        } else if (kmerString.length() - 1 <= MAX_3LONG_ENCODE) {
//            return new PairMer3LongEncoded(leftClip, core, rightClip);
//        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
//            return new PairMer4LongEncoded(leftClip, core, rightClip);
//        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
//            return new PairMer5LongEncoded(leftClip, core, rightClip);
//        } else {
//            return new PairMerIntArrEncoded(leftClip, core, rightClip);
//        }
//    }
}
