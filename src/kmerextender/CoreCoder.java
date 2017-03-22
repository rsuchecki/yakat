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
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.InputReaderProducer;
import shared.Reporter;

/**
 * CAREFUL! - using all bits in signed fields, may cause errors if used for
 * comparisons without decoding (for example lex order of core vs its RC)
 *
 *
 *
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class CoreCoder {

    public static int[] encodeCoreIntArray(String kmerString) {
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

    public static long[] encodeCoreLongArray(String coreString) {
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

    public static void main(String args[]) {
//        String input = "CACCCGGTTAATATAT";        
//        System.err.println(input+" len="+input.length());
//        int encodeCore = encodeCore(input);
//        System.err.println(decodeCore(input.length(), encodeCore));
//        int x =0;
        new CoreCoder();
    }

    public CoreCoder() {

        int prefixLen = 2;
        int prefixMax = (int) Math.pow(4, prefixLen);
        int prefixpLen = 2;
        int prefixpMax = (int) Math.pow(4, prefixpLen);

        PairMer[][] arr = new PairMer[prefixMax][prefixpMax];
        Random random = new Random();
        long count = 0;
        long dups = 0;
        long diffs = 0;
        int randomInts = prefixMax * prefixpMax;

        for (int i = 0; i < randomInts; i++) {
//            CharSequence kmerSequence = decodeCore(16, random.nextInt(1024));
            CharSequence kmerSequence = decodeCore(16, random.nextInt());
            int prefix = encodeCore(kmerSequence.subSequence(0, prefixLen));
            int prefixp = encodeCore(kmerSequence.subSequence(prefixLen, prefixLen + prefixpLen));
            CharSequence suffixString = kmerSequence.subSequence(prefixLen + prefixpLen, kmerSequence.length());
//            int suffix = encodeCore(kmerSequence.subSequence(prefixLen + prefixpLen, kmerSequence.length()));
            //elem not seen before, store
            find:
            if (arr[prefix][prefixp] == null) {
                try {
                    arr[prefix][prefixp] = new PairMerIntArrEncoded(kmerSequence, prefixLen + prefixpLen, kmerSequence.length()-1, false, 1);                    
                    count++;
                } catch (NonACGTException ex) {
                    ex.printStackTrace();
                }
            } else {
                PairMer previous = arr[prefix][prefixp];
                do {
                    if (previous.decodeCore(16-prefixLen-prefixpLen).equals(suffixString)) {
                        //same core - TODO - add PM to core                    
                        dups++;
                        break find;
                    } else {
//                        System.err.println(previous.getSuffix() + " != " + suffix);
                    }

                } while ((previous = previous.getNextPairMer()) != null);
                try {
                    //replace first elem in linked list
                    arr[prefix][prefixp] = new PairMerIntArrEncodedWithRef(kmerSequence, prefixLen + prefixpLen, kmerSequence.length()-1, false, 1, arr[prefix][prefixp]);
                    diffs++;
                } catch (NonACGTException ex) {
                    ex.printStackTrace();
                }
            }
        }
        System.err.println(randomInts + " tried, loaded " + count + ", " + dups + " duplicated and " + diffs + " different at given prefix, total = " + (count + dups + diffs));

        int lev1load = 0;
        long lev2load = 0;
        for (int i = 0; i < arr.length; i++) {
            CharSequence prefixSeq = decodeCore(prefixLen, i);
            if (arr[i] != null) {
                lev1load++;
                int localCount = 0;
                for (int j = 0; j < arr[i].length; j++) {
                    CharSequence prefixpSeq = decodeCore(prefixpLen, j);
                    if (arr[i][j] != null) {
                        lev2load++;
                        localCount++;
                        int local = 0;
                        PairMer p = arr[i][j];
                        do {
                            local++;
                        }while ((p = p.getNextPairMer()) != null);
                        System.err.println(prefixSeq+"_"+prefixpSeq+"_"+arr[i][j].decodeCore(16-prefixLen-prefixpLen)+" "+local);                        
                    }
                }
                System.err.println(prefixSeq+ " total "+localCount);
            }
        }
        System.err.println("load level1 = " + ((double) lev1load) / arr.length + " load level2 = " + (double) lev2load / arr[0].length / arr.length);
    }

//    private class PMO {
//
//        private final int suffix;
//
//        public PMO(int suffix) {
//            this.suffix = suffix;
//        }
//
//        public int getSuffix() {
//            return suffix;
//        }
//
//        public PMO getNextPairMer() {
//            return null;
//        }
//
//    }
//
//    private class PMO1 extends PMO {
//
//        private final PMO next;
//
//        public PMO1(int suffix, PMO next) {
//            super(suffix);
//            this.next = next;
//        }
//
//        @Override
//        public PMO getNextPairMer() {
//            return next;
//        }
//
//    }
}
