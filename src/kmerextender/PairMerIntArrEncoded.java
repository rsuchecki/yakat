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

import java.text.NumberFormat;
import java.util.ArrayList;
import shared.SequenceOps;
import shared.Reporter;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

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
     * @throws kmerextender.NonACGTException
     */
    public PairMerIntArrEncoded(CharSequence sequence, int from, int to, boolean frontClip, int freq) throws NonACGTException {
        addFirstKmer(sequence, from, to, frontClip, freq);
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
     * @throws kmerextender.NonACGTException
     */
    public PairMerIntArrEncoded(CharSequence kmerCoreOnly) throws NonACGTException {
//        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
        encodeCoreCanonical(kmerCoreOnly, 0, kmerCoreOnly.length() - 1);
//        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
    }

    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param encodedKmerCoreOnly
     * @param coreLen
     */
    public PairMerIntArrEncoded(int[] encodedKmerCoreOnly, int coreLen) {
//        String decodeCore = decodeCore(coreLen, encodedKmerCoreOnly);

        storeCanonical(encodedKmerCoreOnly, recodeInReverseComplement(coreLen, encodedKmerCoreOnly, false));
//        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
//        encodeCoreCanonical(kmerCoreOnly, 0, kmerCoreOnly.length() - 1);
//        tmpCore = SequenceOps.getCanonical(kmerCoreOnly);
    }
    
    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param encodedKmerCoreOnly
     * @param coreLen
     */
    public PairMerIntArrEncoded(int[] encodedKmerCoreOnly, int coreLen, boolean canonical) {
//        String decodeCore = decodeCore(coreLen, encodedKmerCoreOnly);
        if(canonical) {
            storeCanonical(encodedKmerCoreOnly, recodeInReverseComplement(coreLen, encodedKmerCoreOnly, false));
        } else {
            kmerCoreBitsArray = encodedKmerCoreOnly;  
        }
//        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
//        encodeCoreCanonical(kmerCoreOnly, 0, kmerCoreOnly.length() - 1);
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
    protected final void addFirstKmer(CharSequence sequence, int from, int to, boolean frontClip, int freq) throws NonACGTException {
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

//            incrementStoredCount();
            incrementStoredCount(hasLeftClip(), freq);
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
//        return decodeCore(44).hashCode();
    }

    @Override
    public int compareTo(PairMerIntArrEncoded anotherKmer) {
//        return decodeCore(44).compareTo(anotherKmer.decodeCore(44));
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

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            sb.append(kmerCoreBitsArray[i]).append(" ");            
        }
        return sb.toString();
    }
    
    /**
     *
     * @param forward
     * @param reverseComplement
     * @return true if stored in forward orientation, false if in
     * reverse-complement
     */
    public final boolean storeCanonical(int[] forward, int[] reverseComplement) { //, String tmpCore) {

//        String canonical = SequenceOps.getCanonical(tmpCore);
        for (int i = 0; i < forward.length; i++) {
            if (forward[i] < reverseComplement[i]) {
                break;
            } else if (forward[i] > reverseComplement[i]) {
                kmerCoreBitsArray = reverseComplement;
//                if(canonical.equals(tmpCore)) {
//                    System.err.println("Cannonical storage error "+tmpCore+" "+canonical);
//                }
                return false;
            }
        }
//        if(!canonical.equals(tmpCore)) {
//            System.err.println("Cannonical storage error "+tmpCore+" "+canonical);
//        }
        kmerCoreBitsArray = forward;
        return true;
    }

    public final boolean isCanonical(int[] encoded, int k) {
        int[] reverseComplement = recodeInReverseComplement(k - 1, encoded, false);
        for (int i = 0; i < encoded.length; i++) {
            if (encoded[i] < reverseComplement[i]) {
                break;
            } else if (encoded[i] > reverseComplement[i]) {
                return false;
            }
        }
        return true;
    }

    public int[] whichCanonical(int[] encoded, int[] reverseComplement, int k) {
        for (int i = 0; i < encoded.length; i++) {
            if (encoded[i] < reverseComplement[i]) {
                break;
            } else if (encoded[i] > reverseComplement[i]) {
                return reverseComplement;
            }
        }
        return encoded;
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

    private int[] recodeInReverseComplement(int encodedSequenceLength, int kmerCoreBitsArray[], boolean debug) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        StringBuilder sb = new StringBuilder();
        int[] rc = new int[kmerCoreBitsArray.length];
        int firstReadChunkLength = encodedSequenceLength * 2 % INT_LENGTH;
        int writeBit = 0;
        int writeChunk = 0; //= rc.length - 1;
        if (debug) {
            System.err.println("Read: ");
        }
        for (int i = kmerCoreBitsArray.length - 1; i > -1; --i) {
            int stopReadingBitsAt = INT_LENGTH;
            int stopWritingBitsAt = INT_LENGTH;
            if (i == 0 && firstReadChunkLength != 0) {
                stopReadingBitsAt = firstReadChunkLength;
            }
            if (writeChunk == 0 && firstReadChunkLength != 0) {
                stopWritingBitsAt = firstReadChunkLength;
            }
            if (debug) {
                System.out.println("Reading from i=" + i + " " + kmerCoreBitsArray[i] + " stop at " + stopReadingBitsAt);
            }
            for (int j = 0; j < stopReadingBitsAt; j += 2) {
                int b2 = (kmerCoreBitsArray[i] & (1 << j)) >> j;
                int b1 = (kmerCoreBitsArray[i] & (1 << j + 1)) >> j + 1;
//                if(debug) {
//                    int x =0;
//                }
//                try {
                rc[writeChunk] <<= 1;
//                } catch (IndexOutOfBoundsException e) {
//                    String forwardCore = decodeCore(encodedSequenceLength, kmerCoreBitsArray);
////                    String reverseCore = decodeCore(encodedSequenceLength, reverseComplement);
//                    System.err.println(forwardCore);
//                    System.err.println("~");
//                    printPaddedEncoded(encodedSequenceLength + 1, kmerCoreBitsArray);
//                    System.err.println("so far:");
//                    printPaddedEncoded(encodedSequenceLength + 1, rc);
//                    
//                    e.printStackTrace();
//                    System.exit(1);
//                }
                if (debug) {
                    System.err.print(b1 + "" + b2 + " ");
                }
                if (b1 == 1 && b2 == 1) {
//                    rc[writeChunk]++;
                    rc[writeChunk] <<= 1;
//                    rc[writeChunk]++;
                } else if (b1 == 0 && b2 == 1) {
                    rc[writeChunk]++;
                    rc[writeChunk] <<= 1;
                } else if (b1 == 1 && b2 == 0) {
                    rc[writeChunk] <<= 1;
                    rc[writeChunk]++;
                } else if (b1 == 0 && b2 == 0) {
                    rc[writeChunk]++;
                    rc[writeChunk] <<= 1;
                    rc[writeChunk]++;
                }
//                System.err.println("Written to "+writeBit+" @ "+writeChunk);
                writeBit++;
                if (debug) {
                    System.out.println("Written to complement of " + b1 + "" + b2 + " to " + writeBit + " in " + writeChunk);
                }
                if ((writeChunk == 0 && writeBit >= stopWritingBitsAt / 2) || (writeChunk > 0 && writeBit >= INT_LENGTH / 2)) {
//                    --writeChunk;
                    writeChunk++;
                    writeBit = 0;
                }
            }
        }
        if (debug) {
            System.err.println();
        }
//         String forwardCore = decodeCore(encodedSequenceLength, kmerCoreBitsArray);
//        if (forwardCore.equals("TCAGTAAGGGTTGCCAAAGGGATG") || forwardCore.equals("CATCCCTTTGGCAACCCTTACTGA")) {
//            printPaddedEncoded(encodedSequenceLength+1, kmerCoreBitsArray);
//            System.err.println(decodeCore(encodedSequenceLength, rc));
//            System.err.println("~~");
//            printPaddedEncoded(encodedSequenceLength+1, rc);
//            System.err.println("");
//            System.err.println("");
//
//            int test[] = new int[]{0b101000000000000011, 0b110000000000000000000000000101};
////            int test[] = new int[]{0b101010101010101010,0b101010101010101010101010101010};
//            printPaddedEncoded(25, test);
//            int[] recodeInReverseComplement = recodeInReverseComplement(24, test,true );
//            printPaddedEncoded(25, recodeInReverseComplement);
//
//            System.exit(1);
//        }

        return rc;
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
     * Encode kmer core in it's two possible representations (forward, rev-comp)
     * but store canonical representation only (decided by lex order)
     *
     * @param sequence
     * @param from inclusive
     * @param to inclusive
     * @return true if stored forward, false if RC
     */
    public final boolean encodeCoreCanonical(CharSequence sequence, int from, int to) throws NonACGTException {
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
                        String message = "Failed ecoding kmerstring to int array. Offending char: " + sequence.charAt(position) + " in " + sequence;
                        throw new NonACGTException(message);
//                        System.err.println();
//                        System.err.println("Offending char: " + sequence.charAt(position));
//                        System.err.println("in " + sequence);
//                        System.err.println("....exiting");
//                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }

        return storeCanonical(kmerCoreBitsArray, kmerCoreBitsArrayRC); //, sequence.subSequence(from , to+1).toString());
    }

    /**
     * UNUSED
     *
     * @param k
     * @return
     */
    @Override
    public PairMer getOtherPairmerCoreLeft(int k) {
        
        
        if (!isInvalid() && getStoredCount() == 2) {

//            String decodedCore = decodeCore(k - 1);
            //Reconstruct k-mer 1
//            String s = getClipLeft() + decodedCore.substring(0, decodedCore.length() - 1) + " <- other core1\n";
//            printPaddedEncoded(k, kmerCoreBitsArray);
            //SHIFT BITS RIGHT TWICE TO DROP 2 RIGHTMOST BITS
            //SET 2 LEFTMOST BITS TO THE VALUE OF CLIP-LEFT
            return new PairMerIntArrEncoded(shiftBitsRightAndAddLeftClip(kmerCoreBitsArray, k - 1), k - 1);
//            int[] leftClipCore = shiftBitsRightAndAddLeftClip(kmerCoreBitsArray, k - 1);
//            int[] rightClipCore = shiftBitsLeftAndAddRightClip(kmerCoreBitsArray, k - 1);
//            String[] s = new String []{decodeCore(k-1, leftClipCore),decodeCore(k-1, rightClipCore)};

//            System.err.println("with leftclip:");
//            printPaddedEncoded(k, leftClipCore);
//            System.err.println("with rightclip:");
//            printPaddedEncoded(k, rightClipCore);
//            System.exit(1);
            //Reconstruct k-mer 2
//            s += decodedCore.substring(1) + getClipRight() + " <- other core2";
            //SHIFT LEFT TWICE ADDING THE VALUE OF CLIP-RIGH
            //MASK 2 LEFTMOST BITS
//            int[] rigthClipCore = shiftBitsLeftAndAddRightClip(kmerCoreBitsArray, k-1);
//            return s;
        }
        return null;
    }

    /**
     * UNUSED
     *
     * @param k
     * @return
     */
    @Override
    public PairMer getOtherPairmerCoreRight(int k) {
        if (!isInvalid() && getStoredCount() == 2) {
            return new PairMerIntArrEncoded(shiftBitsLeftAndAddRightClip(kmerCoreBitsArray, k - 1), k - 1);
        }
        return null;
    }

    /**
     * Method intended to extracting "core" which can then be used as a key for
     * finding neighboring pairMer in the implicit de Brujin graph; For example,
     * given a pairMer A_TTC_G (leftClip_core_rightClip) the method should
     * return the encoded representation of ATTC
     *
     * @param kmerCoreBitsArray
     * @param coreLength
     * @return
     */
    private int[] shiftBitsRightAndAddLeftClip(int[] kmerCoreBitsArray, int coreLength) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int[] shiftedBitsWithLeftClip = new int[kmerCoreBitsArray.length];
        //SHIFT RIGHT ACROSS ARRAY
        int droppedShift = INT_LENGTH - 2; //30bits used in int, -2 bits stored 
        int mask = 0b11; //interested in last two bits only
        int dropped = 0; //nothing dropped-off yet
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            shiftedBitsWithLeftClip[i] = kmerCoreBitsArray[i] >> 2;
            //adding 2 most significant bits dropped at i-1, //no effect at i==0 
            shiftedBitsWithLeftClip[i] |= dropped << droppedShift;
            //storing 2 most significant bits dropped at i (this iteration)
            dropped = kmerCoreBitsArray[i] & mask;
        }

//        printPaddedEncoded(coreLength, shiftedBitsWithLeftClip);
        //SETMOST-SIGNIFICANTBITS TO LEFTCLIP VALUE
        int leftClipShift = coreLength * 2 % INT_LENGTH - 2;
        int leftClipMask = encodeNucleotide(getClipLeft()) << leftClipShift;
        shiftedBitsWithLeftClip[0] |= leftClipMask;

//        printPaddedEncoded(coreLength + 1, shiftedBitsWithLeftClip);
        return shiftedBitsWithLeftClip;
    }

    ///SHIFT LEFT TWICE ADDING THE VALUE OF CLIP-RIGH
    //MASK 2 LEFTMOST BITS
    private int[] shiftBitsLeftAndAddRightClip(int[] kmerCoreBitsArray, int coreLength) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int[] shiftedBitsWithRightClip = new int[kmerCoreBitsArray.length];
        //SHIFT LEFT ACROSS ARRAY
        int clearMask = 0b11000000000000000000000000000000; // dec=-1073741824 //0b11 << INT_LENGTH;
        int matchMask = 0b00110000000000000000000000000000;
        int setRightMostTo = encodeNucleotide(getClipRight());
        for (int i = kmerCoreBitsArray.length - 1; i > -1; --i) {
            //shift bits left 
            shiftedBitsWithRightClip[i] = kmerCoreBitsArray[i] << 2;
//            if (kmerCoreBitsArray[i] == 720978088) {
//                System.err.println("Original:");
//                printPaddedEncoded(coreLength + 1, kmerCoreBitsArray);
//                System.err.println("Shifted:");
//                printPaddedEncoded(coreLength + 1, shiftedBitsWithRightClip);
//            }
            //collect dropped bits 
            int droppedHere = kmerCoreBitsArray[i] & matchMask;
            //clear positions 32,31 which must not be used
            shiftedBitsWithRightClip[i] &= ~clearMask;
//            if (decodeCore(coreLength).equals("CACACATACTCACTTTCATGCACACCGCCACAAC")) {
////                System.err.println("Masked with "+mask);
//                System.err.println("Start here:");
//                printPaddedEncoded(coreLength + 1, kmerCoreBitsArray);                
//                printPaddedEncoded(coreLength + 1, shiftedBitsWithRightClip);                
//                System.err.println("Dropped here "+droppedHere);
//                System.err.println("Dropped shft "+(droppedHere >> INT_LENGTH-2));
//                
//            }

            //set two (empty) least significant bits to the value dropped at previous iteration OR to rightClip at i==0
            shiftedBitsWithRightClip[i] |= setRightMostTo;

            //shift and store for next iteration
            setRightMostTo = droppedHere >> INT_LENGTH - 2;
        }
//        if (decodeCore(coreLength).equals("CACACATACTCACTTTCATGCACACCGCCACAAC")) {
//            System.err.println("rightclip=" + getClipRight());
////            System.err.println(decodeCore(coreLength)+" <===");
//            printPaddedEncoded(coreLength + 1, shiftedBitsWithRightClip);                
//        }
//
        //CLEAR 2 MOST-SIGNIFICANT BITS 
        int leftShift = coreLength * 2 % INT_LENGTH - 2;
        int leftOverflowMask = 0b11 << leftShift + 2;
        shiftedBitsWithRightClip[0] &= (~leftOverflowMask);
//        
//        if (decodeCore(coreLength).equals("CACACATACTCACTTTCATGCACACCGCCACAAC")) {
//            System.err.println("rightclip=" + getClipRight());
//            printPaddedEncoded(coreLength + 1, shiftedBitsWithRightClip);                
//            System.exit(1);
//        }
        return shiftedBitsWithRightClip;
    }

    private int encodeNucleotide(char nucleotide) {
        switch (nucleotide) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                //if C : 01
                return 0b01;
            case 'G':
            case 'g':
                //if G : 10
                return 0b10;
            case 'T':
            case 't':
                //if T : 11
                return 0b11;
        }
        return -1;
    }

//    @Override
    public void printPaddedEncoded(int k, int[] kmerCoreBitsArray) {
        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int incompleteChunkLen = (k - 1) * 2 % INT_LENGTH;
        String decodedCore = decodeCore(k - 1, kmerCoreBitsArray);
        StringBuilder sb = new StringBuilder(String.format("%" + (INT_LENGTH - incompleteChunkLen) + "s", " "));
        for (int i = 0; i < decodedCore.length(); i++) {
            sb.append(" " + decodedCore.charAt(i));
            if (sb.length() % INT_LENGTH == 0) {
                sb.append(" ");
                System.err.print(sb.toString());
                sb = new StringBuilder();
            }
        }
        System.err.println();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            System.err.print(leftPad(Integer.toBinaryString(kmerCoreBitsArray[i]), INT_LENGTH, '0') + " ");
        }
        System.err.println();
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            System.err.printf("%" + INT_LENGTH + "s ", kmerCoreBitsArray[i]);
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
    
//    @Override
//    public ArrayList<PairMer> generateSpacedPairMersForMapSplitting(int k, int chunks) {
//        ArrayList<PairMer> pms = new ArrayList<>(chunks);
//        int chunk = 0b00111111111111111111111111111111 / chunks;        
//        int reminder = 0b00111111111111111111111111111111 % chunks;
//        for (int i = chunk; i <= 0b00111111111111111111111111111111-reminder; i+=chunk) {
//            int[] encoded = new int[kmerCoreBitsArray.length];
//            encoded[0] = i;
////            for (int j = 1; j < encoded.length; j++) {
////                encoded[j] = i;// Integer.MAX_VALUE; //0b00111111111111111111111111111111;
////            }
//            pms.add(new PairMerIntArrEncoded(encoded, k-1, false));
//            System.err.println(NumberFormat.getNumberInstance().format(encoded[0])+" "+NumberFormat.getNumberInstance().format(encoded[1])+" "+NumberFormat.getNumberInstance().format(encoded[2]));
//        }
//        return pms;
//    }

    
    @Override
    public PairMer luckyDipPairMer(int k, PairMer another) {
        int[] encoded = new int[kmerCoreBitsArray.length];
        int[] anotherEnc = ((PairMerIntArrEncoded)another).getEncoded();
//        for (int i = 0; i < encoded.length; i++) {
//            encoded[i] = Math.min(kmerCoreBitsArray[i]+kmerCoreBitsArray[i]/2,r.nextInt(0b00111111111111111111111111111111));
//            encoded[i] = Math.min(kmerCoreBitsArray[i]+kmerCoreBitsArray[i]/2,r.nextInt(0b00111111111111111111111111111111));
//            encoded[i] = ThreadLocalRandom.current().nextInt(0b00111111111111111111111111111111);
//        }        
        encoded[0] = ThreadLocalRandom.current().nextInt(kmerCoreBitsArray[0]+1, anotherEnc[0]-1); //
//        encoded[0] = ThreadLocalRandom.current().nextInt(kmerCoreBitsArray[0]+(anotherEnc[0]-kmerCoreBitsArray[0])/4,anotherEnc[0]-1-(anotherEnc[0]-kmerCoreBitsArray[0])/4);
//        int midPoint = (anotherEnc[0]-kmerCoreBitsArray[0])/2;
//         encoded[0] = midPoint; //ThreadLocalRandom.current().nextInt(midPoint/2,midPoint+midPoint/2);
        return new PairMerIntArrEncoded(encoded, k-1, false);
    }
    
    public int[] getEncoded() {
        return kmerCoreBitsArray;
    }
}
