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
package kmermatch;

import shared.SequenceOps;
import kmerextender.*;
import java.util.Arrays;

/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public final class KmerIntArrIntCounter extends Kmer implements Comparable<KmerIntArrIntCounter> {

    //Object overhead 8 B
    private int[] kmerBitsArray;  //12B + len*4 Bytes 
    private int storedCount;  //4B //TODO 
//    private byte storedCount;  //1B
    //then round to multi of 8

    /**
     * Proper constructor
     *
     * @param kmerString
     * @param count
     */
    public KmerIntArrIntCounter(String kmerString, int count) {
        encodeKmer(SequenceOps.getCanonical(kmerString));
        incrementStoredCount(count);
//        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
    }

    public KmerIntArrIntCounter(String kmerString) {        
        encodeKmer(SequenceOps.getCanonical(kmerString));
        incrementStoredCount(1);
    }
     
    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((KmerIntArrIntCounter) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(this.kmerBitsArray);
        return hash;
    }

    @Override
    public int compareTo(KmerIntArrIntCounter anotherKmer) {
        int[] bitArrayAnother = anotherKmer.getKmerBitsArray();
        for (int i = 0; i < kmerBitsArray.length; i++) {
            try {
                if (kmerBitsArray[i] < bitArrayAnother[i]) {
                    return -1;
                } else if (kmerBitsArray[i] > bitArrayAnother[i]) {
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

    public int[] getKmerBitsArray() {
        return kmerBitsArray;
    }

    @Override
    public String decodeKmer(int coreLength) {
        return decodeKmer(coreLength, kmerBitsArray);
    }

    private String decodeKmer(int encodedSequenceLength, int kmerCoreBitsArray[]) {
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

    public final void encodeKmer(String kmerString) {
//        System.err.println("Encoding seq len=" + kmerString.length());

        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = kmerString.length();
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        kmerBitsArray = new int[intsNeeded];
        char[] kmerCharArray = kmerString.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
                kmerBitsArray[currentInt] <<= 1;
                if (kmerCharArray[position] == 'A' || kmerCharArray[position] == 'a') {
                    //if A : 00
                    kmerBitsArray[currentInt] <<= 1;
                } else if (kmerCharArray[position] == 'C' || kmerCharArray[position] == 'c') {
                    //if C : 01 
                    kmerBitsArray[currentInt] <<= 1;
                    kmerBitsArray[currentInt]++;
                } else if (kmerCharArray[position] == 'G' || kmerCharArray[position] == 'g') {
                    //if G : 10
                    kmerBitsArray[currentInt]++;
                    kmerBitsArray[currentInt] <<= 1;
                } else if (kmerCharArray[position] == 'T' || kmerCharArray[position] == 't') {
                    //if T : 11
                    kmerBitsArray[currentInt]++;
                    kmerBitsArray[currentInt] <<= 1;
                    kmerBitsArray[currentInt]++;
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

     protected synchronized void incrementStoredCount(int count) {
        if (storedCount+count < Integer.MAX_VALUE) {
            storedCount += count;
        } else {
            storedCount = Integer.MAX_VALUE;
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
