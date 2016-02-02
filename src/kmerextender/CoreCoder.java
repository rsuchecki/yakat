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

import java.util.Arrays;

/**
 * CAREFUL!  - using all bits in signed fields, 
 * may cause errors if used for comparisons without decoding 
 * (for example lex order of core vs it's RC)
 * 
 * 
 * 
 * 
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class CoreCoder {

    public static int[] encodeCoreInt(String kmerString) {
        int INT_LENGTH = 32; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = kmerString.length();
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        int kmerCoreBitsArray[] = new int[intsNeeded];
        char[] kmerCharArray = kmerString.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                switch (kmerCharArray[position]) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    default:
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
        return kmerCoreBitsArray;
    }

    public static String decodeCore(int encodedSequenceLength, int kmerCoreBitsArray[]) {
        int INT_LENGTH = 32; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % INT_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = INT_LENGTH - 1;
            if (i == 0 && lastChunkLength != 0) {
                startPrintingBitsFrom = lastChunkLength - 1;
            }
            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
                int b1 = kmerCoreBitsArray[i] >> j & 1;
                int b2 = kmerCoreBitsArray[i] >> j - 1 & 1;
//                int b = (kmerCoreBitsArray[i] & (1 << j)) >> j;
//                int m = (kmerCoreBitsArray[i] & (1 << j - 1)) >> j - 1;
                if (b1 == 0 && b2 == 0) {
                    sb.append("A");
                } else if (b1 == 0 && b2 == 1) {
                    sb.append("C");
                } else if (b1 == 1 && b2 == 0) {
                    sb.append("G");
                } else if (b1 == 1 && b2 == 1) {
                    sb.append("T");
                }
            }
        }
        return sb.toString();
    }

    public static long[] encodeCoreLong(String coreString) {
        int LONG_LENGTH = 64; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = coreString.length();
        int longsNeeded = (int) Math.ceil((double) stringLength * 2 / LONG_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = longsNeeded * LONG_LENGTH - stringLength * 2;
        long kmerCoreBitsArray[] = new long[longsNeeded];
        char[] kmerCharArray = coreString.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < LONG_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                switch (kmerCharArray[position]) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    default:
                        System.err.println("Failed ecoding kmerstring to int array....");
                        System.err.println("Offending char: " + kmerCharArray[position]);
                        System.err.println("in " + coreString);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }
//        String decodeCore = decodeCore(stringLength, kmerCoreBitsArray);
//        if (!decodeCore.equals(coreString)) {
////            System.err.println("Error encoding/decoding " + kmerCoreBits1 + " " + kmerCoreBits2);
//            System.err.println(coreString + " <-core");
//            System.err.println(decodeCore + " <-decoded");
//            decodeCore = decodeCore(stringLength, kmerCoreBitsArray);
//        }
        return kmerCoreBitsArray;
    }

    public static String decodeCore(int encodedSequenceLength, long kmerCoreBitsArray[]) {
        int LONG_LENGTH = 64; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % LONG_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = LONG_LENGTH - 1;
            if (i == 0 && lastChunkLength != 0) {
                startPrintingBitsFrom = lastChunkLength - 1;
            }
            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
                long b1 = kmerCoreBitsArray[i] >> j & 1;
                long b2 = kmerCoreBitsArray[i] >> j - 1 & 1;
                if (b1 == 0 && b2 == 0) {
                    sb.append("A");
                } else if (b1 == 0 && b2 == 1) {
                    sb.append("C");
                } else if (b1 == 1 && b2 == 0) {
                    sb.append("G");
                } else if (b1 == 1 && b2 == 1) {
                    sb.append("T");
                }
            }
        }
        return sb.toString();
    }

    public static int compareCores(int[] core, int[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
                if (core[i] < anotherCore[i]) {
                    return -1;
                } else if (core[i] > anotherCore[i]) {
                    return 1;
                }
            } catch (IndexOutOfBoundsException e) {
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }
    
    public static int compareCores(long[] core, long[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
                if (core[i] < anotherCore[i]) {
                    return -1;
                } else if (core[i] > anotherCore[i]) {
                    return 1;
                }
            } catch (IndexOutOfBoundsException e) {
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }
    
    public static int computeHash(long[] core) {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(core);
        return hash;
    }
    
}
