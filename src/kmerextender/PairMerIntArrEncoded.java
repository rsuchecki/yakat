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
import shared.Reporter;
import java.util.Arrays;

/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerIntArrEncoded extends PairMer implements Comparable<PairMerIntArrEncoded>  {

    //Object overhead 8 B
    private int[] kmerCoreBitsArray;  //12B + len*4 Bytes 
//    private char clipLeft = '#';   //2B                     //TODO: encode to 2-3b if sticking to int array
//    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
//    private byte storedCount;  //1B
//    private boolean invalid;  //1B
//    private boolean visited;
    //then round to multi of 8

    
/**
     * Proper constructor
     *
     * @param leftClip
     * @param core
     * @param rightClip
     */
    public PairMerIntArrEncoded(char leftClip, String core, char rightClip) {
        addFirstKmer(leftClip, core, rightClip);
    }

//    /**
//     * Proper constructor
//     *
//     * @param splitMer
//     */
//    public PairMerIntArrEncoded(SplitMer splitMer) {
//        addFirstKmer(splitMer);
////        tmpCore = splitMer.getCore();
//    }
//    
    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param kmerCoreOnly
     */
    public PairMerIntArrEncoded(String kmerCoreOnly) {
        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
//        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
    }

     /**
     * Add first kmer (encoded as a SplitMer)
     *
     * @param leftClip
     * @param core
     * @param rightClip
     */
    protected final void addFirstKmer(char leftClip, String core, char rightClip) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            encodeCore(core);
            if (leftClip != '#') {
                setClipLeft(leftClip);
            }
            if (rightClip != '#') {
                setClipRight(rightClip);
            }
            incrementStoredCount();
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
        }
    }
    
    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMerIntArrEncoded) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(this.kmerCoreBitsArray);
        return hash;
    }

    @Override
    public int compareTo(PairMerIntArrEncoded anotherKmer) {
        int[] bitArrayAnother = anotherKmer.getKmerCoreBitsArray();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            try {
                if (kmerCoreBitsArray[i] < bitArrayAnother[i]) {
                    return -1;
                } else if (kmerCoreBitsArray[i] > bitArrayAnother[i]) {
                    return 1;
                }
            } catch (IndexOutOfBoundsException e) {
//                System.err.println("IndexOutOfBoundsException trying to compare this:");
//                System.err.println(decodeCore(KMER_LENGTH, kmerCoreBitsArray));
////                System.err.println(bitArray.length + " <- lengths -> " + secondBitArray.length);
//                System.err.println("with another:");
//                System.err.println(anotherKmer.decodeCore(KMER_LENGTH, anotherKmer.getKmerCoreBitsArray()));
////                System.err.println(firstBitArrayAnother.length + " <- lengths -> " + secondBitArrayAnother.length);
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }

    public int[] getKmerCoreBitsArray() {
        return kmerCoreBitsArray;
    }

    @Override
    public String decodeCore(int coreLength) {
        return decodeCore(coreLength, kmerCoreBitsArray);
    }

    private String decodeCore(int encodedSequenceLength, int kmerCoreBitsArray[]) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
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

    public final void encodeCore(String kmerString) {
//        System.err.println("Encoding seq len=" + kmerString.length());

        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = kmerString.length();
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        kmerCoreBitsArray = new int[intsNeeded];
        char[] kmerCharArray = kmerString.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                if (kmerCharArray[position] == 'A' || kmerCharArray[position] == 'a') {
                    //if A : 00
                    kmerCoreBitsArray[currentInt] <<= 1;
                } else if (kmerCharArray[position] == 'C' || kmerCharArray[position] == 'c') {
                    //if C : 01 
                    kmerCoreBitsArray[currentInt] <<= 1;
                    kmerCoreBitsArray[currentInt]++;
                } else if (kmerCharArray[position] == 'G' || kmerCharArray[position] == 'g') {
                    //if G : 10
                    kmerCoreBitsArray[currentInt]++;
                    kmerCoreBitsArray[currentInt] <<= 1;
                } else if (kmerCharArray[position] == 'T' || kmerCharArray[position] == 't') {
                    //if T : 11
                    kmerCoreBitsArray[currentInt]++;
                    kmerCoreBitsArray[currentInt] <<= 1;
                    kmerCoreBitsArray[currentInt]++;
                } else {
                    System.err.println("Failed ecoding kmerstring to int array....");
                    System.err.println("Offending char: " + kmerCharArray[position]);
                    System.err.println("in " + kmerString);
                    System.err.println("....exiting");
                    System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }
    }

//    public int compareToAnotherCore(String anotherKmerCore) {
//        int[] bitArrayAnother = encodeCore(anotherKmerCore);
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
////                System.err.println(firstBitArrayAnother.length + " <- lengths -> " + secondBitArrayAnother.length);
//                System.err.println(e.getMessage());
//            }
//        }
//        return 0;
//    }



}
