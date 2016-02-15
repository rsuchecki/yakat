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
import java.util.stream.IntStream;

/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerIntArrEncoded extends PairMer implements Comparable<PairMerIntArrEncoded> {

    //Object overhead 8 B
    private int[] kmerCoreBitsArray;  //12B + len*4 Bytes 
//    private char clipLeft = '#';   //2B                     //TODO: encode to 2-3b if sticking to int array
//    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
//    private byte storedCount;  //1B
//    private boolean invalid;  //1B
//    private boolean visited;
    //then round to multi of 8

    /**
     * Constructor
     *
     * @param sequence
     * @param from
     * @param to inclusive
     * @param frontClip
     */
    public PairMerIntArrEncoded(CharSequence sequence, int from, int to, boolean frontClip) {
        addFirstKmer(sequence, from, to, frontClip);
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
//        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
        encodeCoreCanonical(kmerCoreOnly, 0, kmerCoreOnly.length() - 1);
//        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
    }

//    /**
//     * Add first kmer (encoded as a SplitMer)
//     *
//     * @param leftClip
//     * @param core
//     * @param rightClip
//     */
//    protected final void addFirstKmer(char leftClip, String core, char rightClip) {
//        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
//            encodeCore(core);
//            if (leftClip != '#') {
//                setClipLeft(leftClip);
//            }
//            if (rightClip != '#') {
//                setClipRight(rightClip);
//            }
//            incrementStoredCount();
//        } else {
//            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
//        }
//    }
    /**
     * Add first mer
     *
     * @param sequence
     * @param from
     * @param to inclusive
     * @param frontClip
     */
//    protected final void addFirstKmer(int leftClipAt, CharSequence sequence, int rightClipAt, int coreStart, int coreEnd) {
    protected final void addFirstKmer(CharSequence sequence, int from, int to, boolean frontClip) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            int coreStart = frontClip ? from + 1 : from;
            int coreEnd = frontClip ? to : to - 1;

            boolean storedInForwardOrient = encodeCoreCanonical(sequence, coreStart, coreEnd);
//            int[] coreInRC = encodeCoreRC(sequence, coreStart, coreEnd);

            if (storedInForwardOrient) {
                if (frontClip) {
                    setClipLeft(sequence.charAt(from));
                } else {
                    setClipRight(sequence.charAt(to));
                }
            } else if (frontClip) {
                setClipRight(SequenceOps.complement(sequence.charAt(from)));
            } else {
                setClipLeft(SequenceOps.complement(sequence.charAt(to)));
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

    /**
     *
     * @param forward
     * @param reverseComplement
     * @return true if stored in forward orientation, false if in reverse-complement
     */
    public boolean storeCanonical(int[] forward, int[] reverseComplement) {
        for (int i = 0; i < forward.length; i++) {
            if (forward[i] < reverseComplement[i]) {
                break;
            } else if (forward[i] > reverseComplement[i]) {
                kmerCoreBitsArray = reverseComplement;
                return false;
            }
        }
        kmerCoreBitsArray = forward;
        return true;
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

//    public final void encodeCore(String kmerString) {
////        System.err.println("Encoding seq len=" + kmerString.length());
//
//        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        int stringLength = kmerString.length();
//        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
////        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
//        int currentInt = 0;
//        int position = 0;
//        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
//        kmerCoreBitsArray = new int[intsNeeded];
//        char[] kmerCharArray = kmerString.toCharArray();
//        while (position < stringLength) {
//            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
//                kmerCoreBitsArray[currentInt] <<= 1;
//                switch (kmerCharArray[position]) {
//                    case 'A':
//                    case 'a':
//                        //if A : 00
//                        kmerCoreBitsArray[currentInt] <<= 1;
//                        break;
//                    case 'C':
//                    case 'c':
//                        //if C : 01
//                        kmerCoreBitsArray[currentInt] <<= 1;
//                        kmerCoreBitsArray[currentInt]++;
//                        break;
//                    case 'G':
//                    case 'g':
//                        //if G : 10
//                        kmerCoreBitsArray[currentInt]++;
//                        kmerCoreBitsArray[currentInt] <<= 1;
//                        break;
//                    case 'T':
//                    case 't':
//                        //if T : 11
//                        kmerCoreBitsArray[currentInt]++;
//                        kmerCoreBitsArray[currentInt] <<= 1;
//                        kmerCoreBitsArray[currentInt]++;
//                        break;
//                    default:
//                        System.err.println("Failed ecoding kmerstring to int array....");
//                        System.err.println("Offending char: " + kmerCharArray[position]);
//                        System.err.println("in " + kmerString);
//                        System.err.println("....exiting");
//                        System.exit(1);
//                }
//                positionInInt += 2;
//                position++;
//            }
//            currentInt++;
//            positionInInt = 0;
//        }
//    }
    /**
     * Encode kmer core in it's two possible representations (forward, rev-comp) but store canonical representation only
     * (decided by lex order)
     *
     * @param sequence
     * @param from inclusive
     * @param to inclusive
     * @return true if stored forward, false if RC
     */
    public final boolean encodeCoreCanonical(CharSequence sequence, int from, int to) {
//        System.err.println("Encoding seq len=" + kmerString.length());

        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        int stringLength = kmerString.length();
        int stringLength = to - from + 1;
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
//        int position = 0;
        int position = from;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        int[] kmerCoreBitsArray = new int[intsNeeded];
        int[] kmerCoreBitsArrayRC = new int[intsNeeded];
        int currentIntRc = intsNeeded - 1;
//        char[] kmerCharArray = kmerString.toCharArray();
        int positionInIntRc = 0;//positionInInt;
        while (position <= to) {
            while ((positionInInt < INT_LENGTH) && (position <= to)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                if (positionInIntRc >= INT_LENGTH) {
                    --currentIntRc;
                    positionInIntRc = 0;
                }
                switch (sequence.charAt(position)) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        //11 in Reverse-Complement                        
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        //10 in Reverse-Complement                        
                        positionInIntRc++;
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        //01 in Reverse-Complement                        
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        positionInIntRc++;

                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        //00 in Reverse-Complement                        
                        positionInIntRc += 2;
                        break;
                    default:
                        System.err.println("Failed ecoding kmerstring to int array....");
                        System.err.println("Offending char: " + sequence.charAt(position));
                        System.err.println("in " + sequence);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }
        return storeCanonical(kmerCoreBitsArray, kmerCoreBitsArrayRC);
    }

    /**
     *
     * @param k
     * @return
     */
    @Override
    public String getOtherCores(int k) {
        if (!isInvalid() && getStoredCount() == 2) {

            String decodedCore = decodeCore(k - 1);

            //Reconstruct k-mer 1
            String s = getClipLeft() + decodedCore.substring(0, decodedCore.length() - 1) + " <- other core1+\n";

            printPaddedEncoded(k, kmerCoreBitsArray);

            //SHIFT BITS RIGHT TWICE TO DROP 2 RIGHTMOST BITS
            //SET 2 LEFTMOST BITS TO THE VALUE OF CLIP-LEFT
            int[] leftClipCore = shiftBitsRightAndAddLeftClip(kmerCoreBitsArray, k-1);
           
            
            
            System.exit(1);

            //Reconstruct k-mer 2
            s += decodedCore.substring(1) + getClipRight() + " <- other core2";

            //SHIFT LEFT TWICE ADDING THE VALUE OF CLIP-RIGH
            //MASK 2 LEFTMOST BITS
//            int[] rigthClipCore = shiftBitsLeftAndAddRightClip(kmerCoreBitsArray, k-1);
            return s;
        }
        return null;
    }

    private int[] shiftBitsRightAndAddLeftClip(int[] kmerCoreBitsArray, int coreLength) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int[] shiftedBitsWithLeftClip = new int[kmerCoreBitsArray.length];        
        //SHIFT RIGHT ACROSS ARRAY
        int droppedShift = INT_LENGTH-2; //30bits used in int, -2 bits stored 
        int mask = 0b11; //interested in last two bits
        int dropped = 0;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            shiftedBitsWithLeftClip[i] = kmerCoreBitsArray[i] >> 2;            
            //adding 2 most significant bits dropped at i-1, //no effect at i==0 
            dropped <<= droppedShift; 
            shiftedBitsWithLeftClip[i] |= dropped; 
            //storing 2 most significant bits dropped at i (this iteration)
            dropped = kmerCoreBitsArray[i] & mask;             
        }            
        
        printPaddedEncoded(coreLength, shiftedBitsWithLeftClip);

        //SETMOST-SIGNIFICANTBITS TO LEFTCLIP VALUE
        int leftClipMask = 0b00;
        switch(getClipLeft()) {
            case 'C':
            case 'c':
                //if C : 01
                leftClipMask = 0b01;
                break;
            case 'G':
            case 'g':
                //if G : 10
                leftClipMask = 0b10;
                break;
            case 'T':
            case 't':
                //if T : 11
                leftClipMask = 0b11;
                break;
        }
        int leftClipShift = coreLength * 2 % INT_LENGTH-2;
        leftClipMask <<= leftClipShift;
        shiftedBitsWithLeftClip[0] |= leftClipMask;
        
        printPaddedEncoded(coreLength+1, shiftedBitsWithLeftClip);
        return shiftedBitsWithLeftClip;
    }
    
    
//    @Override
    public void printPaddedEncoded(int k, int[] kmerCoreBitsArray) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int incompleteChunkLen = (k - 1) * 2 % INT_LENGTH;
        String decodedCore = decodeCore(k - 1, kmerCoreBitsArray);
        StringBuilder sb = new StringBuilder(String.format("%" + (INT_LENGTH - incompleteChunkLen) + "s", " "));
        for (int i = 0; i < decodedCore.length(); i++) {
            StringBuilder append = sb.append(" " + decodedCore.charAt(i));
            if (sb.length() % INT_LENGTH == 0) {
                sb.append(" ");
                System.err.print(sb.toString());
                sb = new StringBuilder();
            }
        }
        System.err.println();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            System.err.print(leftPad(Integer.toBinaryString(kmerCoreBitsArray[i]), 30, '0') + " ");
        }
        System.err.println();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            System.err.printf("%30s ", kmerCoreBitsArray[i]);
        }
        System.err.println();
    }

    private CharSequence leftPad(String bitString, int len, char pad) {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < len - bitString.length(); i++) {
            s.append(pad);
        }
        return s.append(bitString);
    }
}
