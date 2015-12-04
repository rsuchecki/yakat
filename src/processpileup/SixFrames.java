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
 *//*
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
package processpileup;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SixFrames {

//    HashMap<Integer, String> frameToCodonMap;
////    HashMap<Integer, String> reverseFrames;
////    ArrayList<String> forwardFrames;
////    ArrayList<String> reverseFrames;
////    ArrayList<String> forwardFramesTranslated;
////    ArrayList<String> reverseFramesTranslated;
//    CodonDictionary codonTranslationDictionary;
//    Sequence s = new Sequence("", "");
//    boolean translate;
////    int sequenceLength;
//
//    public SixFrames(CodonDictionary codonTranslationDictionary, boolean translate) { //, int sequenceLength) {
//        this.codonTranslationDictionary = codonTranslationDictionary;
//        frameToCodonMap = new HashMap<>(6, 1);
////        forwardFrames = new ArrayList<>(3);
////        reverseFrames = new ArrayList<>(3);
//        this.translate = translate;
////        this.sequenceLength = sequenceLength;
////        forwardFramesTranslated = new ArrayList<>(3);
////        reverseFramesTranslated = new ArrayList<>(3);
//    }
//
//    public void addEmptyCodon(int endPosition, int sequenceLenght) {
//        endPosition++; //switch to 1_indexing
//        frameToCodonMap.put(getForwardFrame(endPosition), "?");
//        frameToCodonMap.put(getReverseFrame(endPosition, sequenceLenght), "?");
//    }
//
//    private Integer getForwardFrame(int endPosition) {
//        if (endPosition % 3 == 0) {
//            return 1;
//        } else if ((endPosition - 1) % 3 == 0) {
//            return 2;
//        } else if ((endPosition - 2) % 3 == 0) {
//            return 3;
//        }
//        int plusOne = endPosition + 1;
//        int plusTwo = endPosition + 2;
//        System.err.println(endPosition + " mod 3 =" + endPosition % 3);
//        System.err.println(plusOne + " mod 3 =" + plusOne % 3);
//        System.err.println(plusTwo + " mod 3 =" + plusTwo % 3);
//        return null;
//    }
//
//    private Integer getReverseFrame(int endPosition, int sequenceLenght) {
//        endPosition = sequenceLenght - endPosition ;
//        if (endPosition % 3 == 0) {
//            return -1;
//        } else if ((endPosition - 1) % 3 == 0) {
//            return -2;
//        } else if ((endPosition - 2) % 3 == 0) {
//            return -3;
//        } else {
//            return null;
//        }
//    }
//
//    public void addCodon(String codon, int endPosition, int sequenceLenght) {
//        endPosition++; //switch to 1_indexing
//        //WE ADD the reverse frames i
//        if (translate) {
//            frameToCodonMap.put(getForwardFrame(endPosition), codonTranslationDictionary.translate(codon));
//            frameToCodonMap.put(getReverseFrame(endPosition, sequenceLenght), codonTranslationDictionary.translate(s.getReverseComplementString(codon.toCharArray())));
//        } else {
//            frameToCodonMap.put(getForwardFrame(endPosition), codon);
//            frameToCodonMap.put(getReverseFrame(endPosition,sequenceLenght), s.getReverseComplementString(codon.toCharArray()));
//        }
//    }
//
////    public void addForward(String codon) {
////        forwardFrames.add(codon);
//////        forwardFramesTranslated.add(codonTranslationDictionary.translate(codon));
////    }
////
////    public void addReverse(String codon) {
////        reverseFrames.add(codon);
////    }
////    public ArrayList<String> getForwardFrames() {
////        return forwardFrames;
////    }
////
////    public ArrayList<String> getReverseFrames() {
////        return reverseFrames;
////    }
//    public HashMap<Integer, String> getFrameToCodonMap() {
//        return frameToCodonMap;
//    }

}
