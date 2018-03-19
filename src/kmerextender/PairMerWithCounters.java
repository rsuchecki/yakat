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

import java.util.ArrayList;
import java.util.Arrays;

/**
 * concept structure to hold up to 2 k-mers paired, by having counters 
 * associated with all possible one-base extensions, one could skip separate k-mer counting
 * and easily exclude k-mers of low frequency and possibly also consider extensions 
 * in somewhat ambiguous cases 
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerWithCounters {// implements Comparable<PairMerWithCounters> {

//    private String pairMerString;  //TODO replace, encoding nucl in int[], use sign to indicate if clip stored at first/last/both - ?
//    private String core; //TODO replace, encoding nucl in int[]
//    private int overhangsAndFlagsIndicatingConflict;
    //Object overhead 8 B
    private int[] kmerCoreBitsArray;  //12B + len*4 Bytes 
    private String tmpCore;
    private char clipLeft = '#';   //2B                     //TODO: encode to 2-3b if sticking to int array
    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
    private byte storedCount;  //1B
    private boolean invalid;  //1B
    private boolean visited;
    //then round to multi of 8

    private byte leftA;
    private byte leftC;
    private byte leftG;
    private byte leftT;
    
    private byte rightA;
    private byte rightC;
    private byte rightG;
    private byte rightT;
    
}

//    /**
//     * Proper constructor
//     *
//     * @param splitMer
//     */
//    public PairMerWithCounters(SplitMer splitMer) {
//        addFirstKmer(splitMer);
////        tmpCore = splitMer.getCore();
//    }
//
//    /**
//     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
//     *
//     * @param kmerCoreOnly
//     */
//    public PairMerWithCounters(String kmerCoreOnly) {
//        kmerCoreBitsArray = encode(SequenceOps.getCanonical(kmerCoreOnly));
////        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
//    }
//
//    @Override
//    public boolean equals(Object anotherKmer) {
//        return compareTo((PairMerWithCounters) anotherKmer) == 0;
//    }
//
//    @Override
//    public int hashCode() {
//        int hash = 3;
//        hash = 79 * hash + Arrays.hashCode(this.kmerCoreBitsArray);
//        return hash;
//    }
//
//    @Override
//    public int compareTo(PairMerWithCounters anotherKmer) {
//        int[] bitArrayAnother = anotherKmer.getKmerCoreBitsArray();
//        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
//            try {
//                if (kmerCoreBitsArray[i] < bitArrayAnother[i]) {
//                    return -1;
//                } else if (kmerCoreBitsArray[i] > bitArrayAnother[i]) {
//                    return 1;
//                }
//            } catch (IndexOutOfBoundsException e) {
////                System.err.println("IndexOutOfBoundsException trying to compare this:");
////                System.err.println(decodeCore(KMER_LENGTH, kmerCoreBitsArray));
//////                System.err.println(bitArray.length + " <- lengths -> " + secondBitArray.length);
////                System.err.println("with another:");
////                System.err.println(anotherKmer.decodeCore(KMER_LENGTH, anotherKmer.getKmerCoreBitsArray()));
//////                System.err.println(firstBitArrayAnother.length + " <- lengths -> " + secondBitArrayAnother.length);
//                System.err.println(e.getMessage());
//
//            }
//        }
//        return 0;
//    }
////    public boolean addKmer(String kmerString, boolean frontClip, int overlapLength) {
//
//    /**
//     * Add first kmer (encoded as a SplitMer)
//     *
//     * @param split : SplitMer representation of a k-mer
//     */
//    private void addFirstKmer(SplitMer split) {
////        if (split.getCore().equals("AAAACGTAAGAGTATCATTCAGCTGTTATG")) {
////            int x = 0;
////        }
////        history.add("[1]\t" + split.getLeftClip() + "_" + split.getCore() + "_" + split.getRightClip());
////        if (invalid || storedCount > 1) {
////            invalid = true;
////            return false;
////        }
////        SplitMer split = new SplitMer(kmerString, frontClip, overlapLength);
//        if (storedCount == 0) {        //If this is the first of the two k-mers that could be stored
//            kmerCoreBitsArray = encode(split.getCore());
//            if (!split.getLeftClip().isEmpty()) {
//                clipLeft = split.getLeftClipChar();
//            }
//            if (!split.getRightClip().isEmpty()) {
//                clipRight = split.getRightClipChar();
//            }
//            incrementStoredCount();
//        } else {
//            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!");
//        }
////            else if (storedCount == 1) {                        //there already is one k-mer stored here
////            if ((clipLeft == '#' && split.getLeftClip().isEmpty()) || (clipRight == '#' && split.getRightClip().isEmpty())) {
////                invalid = true;
////            } else if (clipLeft == '#' && !split.getLeftClip().isEmpty()) {
////                clipLeft = split.getLeftClipChar();
////            } else if (clipRight == '#' && !split.getRightClip().isEmpty()) {
////                clipRight = split.getLeftClipChar();
////            } else {
////                Reporter.report("[ERROR]", "Unexpected error [2], addKmer() at " + this.getClass().getSimpleName());
////            }
////        } else { //already 2 or more k-mers stored, so extension not possible!
////            invalid = true;
////        }
//    }
//
//    /**
//     * It is assumed that the input PairMer matches this one (another.core == this.core)
//     * @param another : newly generated PairMer holding a single k-mer
//     * @param inputKmersUnique : set true if no duplicate k-mers expected
//     */
//    public synchronized void addKmerSynchronized(PairMerWithCounters another, boolean inputKmersUnique) {
//        if (invalid || (inputKmersUnique && storedCount > 1) || (!inputKmersUnique && storedCount > 2)) { //if already invalid PairMer  or more than second kmer being added
//            invalid = true;
//        } else {
//            //If input k-mers are non-unique, ie, a k-mer may appear more than once, we need to ensure that we ignore it            
//            if (!inputKmersUnique) {
//                if ((hasLeftClip() && clipLeft == another.getClipLeft()) || (hasRightClip() && clipRight == another.getClipRight())) {
//                    //fine, new k-mer is identical to a k-mer already stored
////                    System.err.println("Adding same kmer again");
//                } else { //it is a different k-mer
//                    if ((hasLeftClip() && another.hasLeftClip()) || (hasRightClip() && another.hasRightClip())) {
//                        invalid = true;
//                    } else if (!hasLeftClip() && another.hasLeftClip()) {
//                        clipLeft = another.getClipLeft();
//                    } else if (!hasRightClip() && another.hasRightClip()) {
//                        clipRight = another.getClipRight();
//                    } else {
//                        Reporter.report("[BUG?]", "Unexpected [2], addKmer() at " + this.getClass().getSimpleName());
//                    }
//                    incrementStoredCount();
//                }
//            } else { //i.e. input k-mers are unique (no duplicates which we need to ignore)
//                if (hasBothClips()) {
//                    invalid = true;
//                } else if((hasLeftClip() && another.hasLeftClip()) || (hasRightClip() && another.hasRightClip())) {
//                    invalid = true;
//                } else if (!hasLeftClip() && another.hasLeftClip()) {
//                    clipLeft = another.getClipLeft();
//                } else if (!hasRightClip() && another.hasRightClip()) {
//                    clipRight = another.getClipRight();
//                } else {
//                    Reporter.report("[BUG?]", "Unexpected [3], addKmer() at " + this.getClass().getSimpleName());
//                }
//                incrementStoredCount();
//            }
//        }
//    }
//
//    public boolean hasLeftClip() {
//        return clipLeft != '#';
//    }
//
//    public boolean hasRightClip() {
//        return clipRight != '#';
//    }
//
//    public boolean hasBothClips() {
//        return hasLeftClip() && hasRightClip();
//    }
//
//    private void incrementStoredCount() {
//        if (storedCount < Byte.MAX_VALUE) {
//            storedCount++;
//        }
//    }
//
//    public int[] getKmerCoreBitsArray() {
//        return kmerCoreBitsArray;
//    }
//
//    public char getClipLeft() {
//        return clipLeft;
//    }
//
//    public char getClipRight() {
//        return clipRight;
//    }
//
//    public String decodeCore(int coreLength) {
//        return decodeCore(coreLength, kmerCoreBitsArray);
//    }
//
//    private String decodeCore(int encodedSequenceLength, int kmerCoreBitsArray[]) {
//        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        StringBuilder sb = new StringBuilder();
//        int lastChunkLength = encodedSequenceLength * 2 % INT_LENGTH;
//        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
//            int startPrintingBitsFrom = INT_LENGTH - 1;
//            if (i == 0 && lastChunkLength != 0) {
//                startPrintingBitsFrom = lastChunkLength - 1;
////            } else {
////                System.out.println("Printing bits from "+startPrintingBitsFrom);
//            }
//            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
//                int b = (kmerCoreBitsArray[i] & (1 << j)) >> j;
//                int m = (kmerCoreBitsArray[i] & (1 << j - 1)) >> j - 1;
//                if (b == 0 && m == 0) {
//                    sb.append("A");
//                } else if (b == 0 && m == 1) {
//                    sb.append("C");
//                } else if (b == 1 && m == 0) {
//                    sb.append("G");
//                } else if (b == 1 && m == 1) {
//                    sb.append("T");
//                }
//            }
//        }
//        return sb.toString();
//    }
//
//    private int[] encode(String kmerString) {
////        System.err.println("Encoding seq len=" + kmerString.length());
//
//        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        int stringLength = kmerString.length();
//        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
////        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
//        int currentInt = 0;
//        int position = 0;
//        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
//        int[] bitsArray = new int[intsNeeded];
//        char[] kmerCharArray = kmerString.toCharArray();
//        while (position < stringLength) {
//            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
//                bitsArray[currentInt] <<= 1;
//                if (kmerCharArray[position] == 'A' || kmerCharArray[position] == 'a') {
//                    //if A : 00
//                    bitsArray[currentInt] <<= 1;
//                } else if (kmerCharArray[position] == 'C' || kmerCharArray[position] == 'c') {
//                    //if C : 01 
//                    bitsArray[currentInt] <<= 1;
//                    bitsArray[currentInt]++;
//                } else if (kmerCharArray[position] == 'G' || kmerCharArray[position] == 'g') {
//                    //if G : 10
//                    bitsArray[currentInt]++;
//                    bitsArray[currentInt] <<= 1;
//                } else if (kmerCharArray[position] == 'T' || kmerCharArray[position] == 't') {
//                    //if T : 11
//                    bitsArray[currentInt]++;
//                    bitsArray[currentInt] <<= 1;
//                    bitsArray[currentInt]++;
//                } else {
//                    System.err.println("Failed ecoding kmerstring to int array....");
//                    System.err.println("Offending char: " + kmerCharArray[position]);
//                    System.err.println("in " + kmerString);
//                    System.err.println("....exiting");
//                    System.exit(1);
//                }
//                positionInInt += 2;
//                position++;
//            }
//            currentInt++;
//            positionInInt = 0;
//        }
//        return bitsArray;
//    }
//
//    public int compareToAnotherCore(String anotherKmerCore) {
//        int[] bitArrayAnother = encode(anotherKmerCore);
//        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
//            try {
//                if (kmerCoreBitsArray[i] < bitArrayAnother[i]) {
//                    return -1;
//                } else if (kmerCoreBitsArray[i] > bitArrayAnother[i]) {
//                    return 1;
//                }
//            } catch (IndexOutOfBoundsException e) {
////                System.err.println("IndexOutOfBoundsException trying to compare this:");
////                System.err.println(decodeCore(KMER_LENGTH, kmerCoreBitsArray));
//////                System.err.println(bitArray.length + " <- lengths -> " + secondBitArray.length);
////                System.err.println("with another:");
////                System.err.println(decodeCore(KMER_LENGTH, bitArrayAnother));
//////                System.err.println(firstBitArrayAnother.length + " <- lengths -> " + secondBitArrayAnother.length);
//                System.err.println(e.getMessage());
//            }
//        }
//        return 0;
//    }
//
//    public boolean isInvalid() {
//        return invalid;
//    }
//
//    public byte getStoredCount() {
//        return storedCount;
//    }
//
//    public boolean isVisited() {
//        return visited;
//    }
//
//    public void setVisited() {
//        this.visited = true;
//    }
//
//    public String getPairMerString(int k) {
//        StringBuilder sb = new StringBuilder();
//        sb.append(getClipLeft());
//        sb.append(decodeCore(k - 1));
//        sb.append(getClipRight());
//        return sb.toString();
//    }
//
//}
