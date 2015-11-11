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
 * One of the two representations of a k-mer, with one of the ends (typically 1 base)
 * stored separately from the reminder (core - typically k-1 bases). 
 * The object has the core field filled and either the leftCLip or the rightCLip
 * Only used as a wrapper class when generating more compact PairMer 
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SplitMer {

        private final String leftClip;
        private final String core;
        private final String rightClip;

        public SplitMer(String kmerString, boolean frontClip, int overlapLength) {
            //SPLIT THE INPUT INTO CORE AND CLIP
            String coreTmp;
            String clip;
            int len = kmerString.length();
            if (frontClip) {
                coreTmp = kmerString.substring(len - overlapLength);
                clip = kmerString.substring(0, len - overlapLength);
            } else {
                coreTmp = kmerString.substring(0, overlapLength);
                clip = kmerString.substring(overlapLength);
            }

            //ORIENTATE CORE AND CLIP BASED ON LEX ORDER OF CORE AND ITS REV-COMP
            String coreRC = SequenceOps.getReverseComplementString(coreTmp);
            if (coreRC.compareTo(coreTmp) < 0) {
                //REV_COMP = TRUE
                this.core = coreRC;
                if (frontClip) {
                    this.leftClip = "";
                    this.rightClip = SequenceOps.getReverseComplementString(clip);
                } else {
                    this.leftClip = SequenceOps.getReverseComplementString(clip);
                    this.rightClip = "";
                }
            } else {
                this.core = coreTmp;
                if (frontClip) {
                    this.leftClip = clip;
                    this.rightClip = "";
                } else {
                    this.leftClip = "";
                    this.rightClip = clip;
                }
            }
        }

        public String getLeftClip() {
            return leftClip;
        }

        public char getLeftClipChar() {
            return leftClip.charAt(0);
        }

        public String getCore() {
            return core;
        }

        public String getRightClip() {
            return rightClip;
        }

        public char getRightClipChar() {
            return rightClip.charAt(0);
        }
    }