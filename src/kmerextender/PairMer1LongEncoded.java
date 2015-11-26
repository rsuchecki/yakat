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

/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMer1LongEncoded extends PairMer implements Comparable<PairMer1LongEncoded> {

    private long kmerCoreBits;

    /**
     * Proper constructor
     *
     * @param splitMer
     */
    public PairMer1LongEncoded(SplitMer splitMer) {
        addFirstKmer(splitMer);
    }

    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param kmerCoreOnly
     */
    public PairMer1LongEncoded(String kmerCoreOnly) {
        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
    }

    /**
     * Add first kmer (encoded as a SplitMer)
     *
     * @param split : SplitMer representation of a k-mer
     */
    private void addFirstKmer(SplitMer split) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored

            encodeCore(split.getCore());
            if (!split.getLeftClip().isEmpty()) {
                setClipLeft(split.getLeftClipChar());
            }
            if (!split.getRightClip().isEmpty()) {
                setClipRight(split.getRightClipChar());
            }
            incrementStoredCount();
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName() );
        }
    }

    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMer1LongEncoded) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
//        int hash = 3;
//        hash = 79 * hash + Long.valueOf(kmerCoreBits).hashCode();
//        return hash;
        return (int) (kmerCoreBits ^ (kmerCoreBits >>> 32)); //as in java Long hashCode()
    }

    @Override
    public int compareTo(PairMer1LongEncoded anotherKmer) {
        long bitsAnother = anotherKmer.getKmerCoreBits();
        if (kmerCoreBits < bitsAnother) {
            return -1;
        } else if (kmerCoreBits > bitsAnother) {
            return 1;
        }
        return 0;
    }

    public long getKmerCoreBits() {
        return kmerCoreBits;
    }

    @Override
    public String decodeCore(int coreLength) {
        return decodeCore(coreLength, kmerCoreBits);
    }

    private String decodeCore(int encodedSequenceLength, long kmerCoreBits) {
//        int LONG_LENGTH = 62; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
//        int lastChunkLength = encodedSequenceLength * 2 % LONG_LENGTH;
//        int startPrintingBitsFrom = LONG_LENGTH - 1;
//        if ( lastChunkLength != 0) {
        int startPrintingBitsFrom = encodedSequenceLength * 2 - 1;
//            } else {
//                System.out.println("Printing bits from "+startPrintingBitsFrom);
//        }
//        if (kmerCoreBits == 170160154469L) {
//            String toBinaryString = Long.toBinaryString(kmerCoreBits);
//            System.err.println(toBinaryString);
//            for (int i = 64; i > -1; i--) {
//                System.err.printf("%2d %d\n",i,(kmerCoreBits >>> i & 1));
//            }
//            System.err.println("DONE");
//        }

        for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
            long twoBits = kmerCoreBits >> j-1 & 3;
            if(twoBits == 0) {
                sb.append("A");                
            } else if(twoBits == 1) {
                sb.append("C");                
            } else if(twoBits == 2) {
                sb.append("G");                
            } else if(twoBits == 3) {
                sb.append("T");                
            }
//            long b1 = kmerCoreBits >> j & 1;
//            long b2 = kmerCoreBits >> j - 1 & 1;
//            if (b1 == 0 && b2 == 0) {
//                sb.append("A");
//            } else if (b1 == 0 && b2 == 1) {
//                sb.append("C");
//            } else if (b1 == 1 && b2 == 0) {
//                sb.append("G");
//            } else if (b1 == 1 && b2 == 1) {
//                sb.append("T");
//            }
        }
        return sb.toString();
    }

    public final void encodeCore(String coreString) {
//        System.err.println("Encoding seq len=" + kmerString.length());
//        int LONG_LENGTH = 62; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = coreString.length();
//        int currentInt = 0;
        int position = 0;
//        kmerCoreBits = new int[intsNeeded];
        char[] kmerCharArray = coreString.toCharArray();
        while (position < stringLength) {
            kmerCoreBits <<= 1;
            if (kmerCharArray[position] == 'A' || kmerCharArray[position] == 'a') {
                //if A : 00
                kmerCoreBits <<= 1;
            } else if (kmerCharArray[position] == 'C' || kmerCharArray[position] == 'c') {
                //if C : 01 
                kmerCoreBits <<= 1;
                kmerCoreBits++;
            } else if (kmerCharArray[position] == 'G' || kmerCharArray[position] == 'g') {
                //if G : 10
                kmerCoreBits++;
                kmerCoreBits <<= 1;
            } else if (kmerCharArray[position] == 'T' || kmerCharArray[position] == 't') {
                //if T : 11
                kmerCoreBits++;
                kmerCoreBits <<= 1;
                kmerCoreBits++;
            } else {
                System.err.println("Failed ecoding kmerstring to long....");
                System.err.println("Offending char: " + kmerCharArray[position]);
                System.err.println("in " + coreString);
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
    }
}
