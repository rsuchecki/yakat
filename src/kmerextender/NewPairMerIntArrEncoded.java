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

import shared.Reporter;

/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class NewPairMerIntArrEncoded extends PairMer implements Comparable<NewPairMerIntArrEncoded>{

    //Object overhead 8 B
    private int[] kmerCoreReminderBitsArray;  //12B + len*4 Bytes 
//    private char clipLeft = '#';   //2B                     //TODO: encode to 2-3b if sticking to int array
//    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
//    private byte storedCount;  //1B
//    private boolean invalid;  //1B
//    private boolean visited;
    //then round to multi of 8

    /**
     * Constructor
     *
     * @param canonicalCoreReminder
     * @param clip
     * @param frontClip
     * @param freq
     * @throws kmerextender.NonACGTException
     */
    public NewPairMerIntArrEncoded(CharSequence canonicalCoreReminder, char clip, boolean frontClip, int freq) throws NonACGTException {
        addFirstKmer(canonicalCoreReminder, clip, frontClip, freq);
    }
    
    public NewPairMerIntArrEncoded(int[] encodedCanonicalCoreReminded, char clip, boolean frontClip, int freq)  {
        addFirstKmer(encodedCanonicalCoreReminded, clip, frontClip, freq);
    }

//    /**
    /**
     * Add first mer
     *
     * @param canonicalCoreReminder
     * @param clip
     * @param frontClip
     * @param freq
     * @throws kmerextender.NonACGTException
     */
//    protected final void addFirstKmer(int leftClipAt, CharSequence sequence, int rightClipAt, int coreStart, int coreEnd) {
    protected final void addFirstKmer(CharSequence canonicalCoreReminder, char clip, boolean frontClip, int freq) throws NonACGTException {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            if (frontClip) {
                setClipLeft(clip);
            } else {
                setClipRight(clip);
            }
            kmerCoreReminderBitsArray = CoreCoder.encodeCoreIntArray(canonicalCoreReminder);
            incrementStoredCount(hasLeftClip(), freq);
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
        }
    }

    /**
     *
     * @param encodedCanonicalCoreReminded
     * @param clip
     * @param frontClip
     * @param freq
     */
    protected final void addFirstKmer(int[] encodedCanonicalCoreReminded, char clip, boolean frontClip, int freq) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            if (frontClip) {
                setClipLeft(clip);
            } else {
                setClipRight(clip);
            }
            kmerCoreReminderBitsArray = encodedCanonicalCoreReminded;
            incrementStoredCount(hasLeftClip(), freq);
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
        }
    }

    @Override
    public String decodeCore(int coreLength) {
        return decodeCore(coreLength, kmerCoreReminderBitsArray);
    }

    private String decodeCore(int encodedSequenceLength, int kmerCoreBitsArray[]) {
        int INT_LENGTH = 32; //using sign bit as well
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % INT_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = INT_LENGTH - 1;
            if (i == 0 && lastChunkLength != 0) {
                startPrintingBitsFrom = lastChunkLength - 1;
//            } else {
//                System.out.println("Printing bits from "+startPrintingBitsFrom);
            }
            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
                int b = (kmerCoreBitsArray[i] & (1 << j)) >> j;
                int m = (kmerCoreBitsArray[i] & (1 << j - 1)) >> j - 1;
                if (b == 0 && m == 0) {
                    sb.append("A");
                } else if (b == 0 && m == 1) {
                    sb.append("C");
                } else if (b == 1 && m == 0) {
                    sb.append("G");
                } else if (b == 1 && m == 1) {
                    sb.append("T");
                }
            }
        }
        return sb.toString();
    }

    @Override
    public int compareTo(NewPairMerIntArrEncoded anotherMer) {
        return CoreCoder.compareCores(getBitFields(), anotherMer.getBitFields());
    }

    public int[] getBitFields() {
        return kmerCoreReminderBitsArray;
    }

}
