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
package kextender;


import java.util.Arrays;

/**
 *
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class CoreCoder {

    public static byte[] encodeCoreByteArray(CharSequence kmerSequence) {
        int BYTE_LENGTH = 8; 
        int stringLength = kmerSequence.length();
        int bytesNeeded = (int) Math.ceil((double) stringLength * 2 / BYTE_LENGTH); 
        int currentByte = 0;
        int position = 0;
        int positionInByte = bytesNeeded * BYTE_LENGTH - stringLength * 2;
        byte kmerCoreBitsArray[] = new byte[bytesNeeded];
//        char[] kmerCharArray = kmerSequence.toCharArray();
        while (position < stringLength) {
            while ((positionInByte < BYTE_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentByte] <<= 1;
                switch (kmerSequence.charAt(position)) {
                    case 'A':
                    case 'a':                        
                        kmerCoreBitsArray[currentByte] <<= 1; //if A : 00
                        break;
                    case 'C':
                    case 'c':                        
                        kmerCoreBitsArray[currentByte] <<= 1; //if C : 01
                        kmerCoreBitsArray[currentByte]++;
                        break;
                    case 'G':
                    case 'g':                        
                        kmerCoreBitsArray[currentByte]++; //if G : 10
                        kmerCoreBitsArray[currentByte] <<= 1;
                        break;
                    case 'T':
                    case 't':                        
                        kmerCoreBitsArray[currentByte]++; //if T : 11
                        kmerCoreBitsArray[currentByte] <<= 1;
                        kmerCoreBitsArray[currentByte]++;
                        break;
                    default:
                        System.err.println("Failed encoding kmerstring to byte array....");
                        System.err.println("Offending char: \"" + kmerSequence.charAt(position)+"\"");
                        System.err.println("at position " + position);
                        System.err.println("in " + kmerSequence);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInByte += 2;
                position++;
            }
            currentByte++;
            positionInByte = 0;
        }
        return kmerCoreBitsArray;
    }
    
    public static String decodeCore(int encodedSequenceLength, byte kmerCoreBitsArray[]) {
        int BYTE_LENGTH = 8; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % BYTE_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = BYTE_LENGTH - 1;
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
    
    public static int[] encodeCoreIntArray(CharSequence kmerSequence) {
        int INT_LENGTH = 32; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = kmerSequence.length();
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        int kmerCoreBitsArray[] = new int[intsNeeded];
//        char[] kmerCharArray = kmerSequence.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
//                switch (kmerCharArray[position]) {
                switch (kmerSequence.charAt(position)) {
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
                        System.err.println("Failed ecoding kmerstring to tmp int array....");
                        System.err.println("Offending char: " + kmerSequence.charAt(position));
                        System.err.println("in " + kmerSequence);
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

    public static long[] encodeCoreLongArray(CharSequence coreSequence) {
        int LONG_LENGTH = 64; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = coreSequence.length();
        int longsNeeded = (int) Math.ceil((double) stringLength * 2 / LONG_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = longsNeeded * LONG_LENGTH - stringLength * 2;
        long kmerCoreBitsArray[] = new long[longsNeeded];
//        char[] kmerCharArray = coreSeq.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < LONG_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                switch (coreSequence.charAt(position)) {
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
                        System.err.println("Failed ecoding kmerstring to tmp long array....");
                        System.err.println("Offending char: " + coreSequence.charAt(position));
                        System.err.println("in " + coreSequence);
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

//    public static int compareCores(int[] core, int[] anotherCore) {
//        for (int i = 0; i < core.length; i++) {
//            try {
//                if (core[i] < anotherCore[i]) {
//                    return -1;
//                } else if (core[i] > anotherCore[i]) {
//                    return 1;
//                }
//            } catch (IndexOutOfBoundsException e) {
//                System.err.println(e.getMessage());
//            }
//        }
//        return 0;
//    }
    public static int compareCores(byte[] core, byte[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
                //unisgned comparison to use all bits
                if (core[i] == anotherCore[i]) {
                    //keep going
                } else if ((core[i] < anotherCore[i]) ^ (core[i] < 0) ^ (anotherCore[i] < 0)) {
                    return 1;
                } else {
                    return -1;
                }
            } catch (IndexOutOfBoundsException e) {
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }
    public static int compareCores(int[] core, int[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
                //unisgned comparison to use all bits
                if (core[i] == anotherCore[i]) {
                    //keep going
                } else if ((core[i] < anotherCore[i]) ^ (core[i] < 0) ^ (anotherCore[i] < 0)) {
                    return 1;
                } else {
                    return -1;
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
//                if (core[i] < anotherCore[i]) {
//                (longA < longB) ^ (longA < 0) ^ (longB< 0) ? 1 : -1;
                //unisgned comparison to use all bits
                if (core[i] == anotherCore[i]) {
                    //keep going
                } else if ((core[i] < anotherCore[i]) ^ (core[i] < 0) ^ (anotherCore[i] < 0)) {
                    return 1;
                } else {
                    return -1;
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

    public static CharSequence decodeCore(int encodedSequenceLength, int ecodedCore) {
        StringBuilder sb = new StringBuilder();
        int startPrintingBitsFrom = encodedSequenceLength * 2 - 1;
//        if (kmerCoreBits == 170160154469L) {
//            String toBinaryString = Long.toBinaryString(kmerCoreBits);
//            System.err.println(toBinaryString);
//            for (int i = 64; i > -1; i--) {
//                System.err.printf("%2d %d\n",i,(kmerCoreBits >>> i & 1));
//            }
//            System.err.println("DONE");
//        }

        for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
            int twoBits = ecodedCore >> j - 1 & 3;
            switch (twoBits) {
                case 0:
                    sb.append("A");
                    break;
                case 1:
                    sb.append("C");
                    break;
                case 2:
                    sb.append("G");
                    break;
                case 3:
                    sb.append("T");
                    break;
                default:
                    break;
            }
        }
        return sb;
    }

    public static int encodeCore(CharSequence sequence) {
//        System.err.println("Encoding seq len=" + kmerString.length());
        int stringLength = sequence.length();
//        int currentInt = 0;
        int position = 0;
        int kmerCoreBits = 0;
        while (position < stringLength) {
            kmerCoreBits <<= 1;
            switch (sequence.charAt(position)) {
                case 'A':
                case 'a':
                    //if A : 00
                    kmerCoreBits <<= 1;
                    break;
                case 'C':
                case 'c':
                    //if C : 01
                    kmerCoreBits <<= 1;
                    kmerCoreBits++;
                    break;
                case 'G':
                case 'g':
                    //if G : 10
                    kmerCoreBits++;
                    kmerCoreBits <<= 1;
                    break;
                case 'T':
                case 't':
                    //if T : 11
                    kmerCoreBits++;
                    kmerCoreBits <<= 1;
                    kmerCoreBits++;
                    break;
                default:
                    System.err.println("Failed ecoding kmerstring to long....");
                    System.err.println("Offending char: " + sequence.charAt(position));
                    System.err.println("in " + sequence);
                    System.err.println("....exiting");
                    System.exit(1);
            }
            position++;
        }
//        String decodeCore = decodeCore(stringLength);
//        if (!decodeCore.equals(coreString)) {
//            System.err.println("Error encoding/decoding " + kmerCoreBits);
//            System.err.println(coreString + " <-core");
//            System.err.println(decodeCore + " <-decoded");
//            decodeCore = decodeCore(stringLength);
//        }
        return kmerCoreBits;
    }
}
