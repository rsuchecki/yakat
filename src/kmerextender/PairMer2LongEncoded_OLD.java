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
public class PairMer2LongEncoded_OLD extends PairMer implements Comparable<PairMer2LongEncoded_OLD> {

    private long kmerCoreBits1;
    private long kmerCoreBits2;

    /**
     * Proper constructor
     *
     * @param leftClip
     * @param core
     * @param rightClip
     */
    public PairMer2LongEncoded_OLD(char leftClip, String core, char rightClip, int freq) {
        addFirstKmer(leftClip, core, rightClip, freq);
    }

    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param kmerCoreOnly
     */
    public PairMer2LongEncoded_OLD(String kmerCoreOnly) {
        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
    }

    public final void addFirstKmer(char leftClip, String core, char rightClip, int freq) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            encodeCore(core);
            if (leftClip != '#') {
                setClipLeft(leftClip);
            }
            if (rightClip != '#') {
                setClipRight(rightClip);
            }
            incrementStoredCount(hasLeftClip(), freq);
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
        }
    }

    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMer2LongEncoded_OLD) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        return (int) (kmerCoreBits1 ^ kmerCoreBits2);
    }

    @Override
    public int compareTo(PairMer2LongEncoded_OLD anotherKmer) {
        long bitsAnother = anotherKmer.getKmerCoreBits1();
        if (kmerCoreBits1 < bitsAnother) {
            return -1;
        } else if (kmerCoreBits1 > bitsAnother) {
            return 1;
        } else {
            bitsAnother = anotherKmer.getKmerCoreBits2();
            if (kmerCoreBits2 < bitsAnother) {
                return -1;
            } else if (kmerCoreBits2 > bitsAnother) {
                return 1;
            }
            return 0;
        }
    }

    public long getKmerCoreBits1() {
        return kmerCoreBits1;
    }

    public long getKmerCoreBits2() {
        return kmerCoreBits2;
    }

    @Override
    public String decodeCore(int coreLength) {
        return decodeCore(coreLength, kmerCoreBits1, kmerCoreBits2);
    }

    private String decodeCore(int encodedSequenceLength, long kmerCoreBits1, long kmerCoreBits2) {
        StringBuilder sb = new StringBuilder();
        int startPrintingBitsFrom = encodedSequenceLength - 1;
        for (int j = startPrintingBitsFrom; j > -1; j--) {
            long b1 = kmerCoreBits1 >> j & 1;
            long b2 = kmerCoreBits2 >> j & 1;
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
        return sb.toString();
    }

    public final void encodeCore(String coreString) {
        int stringLength = coreString.length();
        int position = 0;
        char[] kmerCharArray = coreString.toCharArray();
        while (position < stringLength) {
            kmerCoreBits1 <<= 1;
            kmerCoreBits2 <<= 1;
            switch (kmerCharArray[position]) {
            //if A : 00
                case 'A':
                case 'a':
                    break;
                case 'C':
                case 'c':
                    //if C : 01
                    kmerCoreBits2++;
                    break;
                case 'G':
                case 'g':
                    //if G : 10
                    kmerCoreBits1++;
                    break;
                case 'T':
                case 't':
                    //if T : 11
                    kmerCoreBits1++;
                    kmerCoreBits2++;
                    break;
                default:
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
//            System.err.println("Error encoding/decoding " + kmerCoreBits1 + " " + kmerCoreBits2);
//            System.err.println(coreString + " <-core");
//            System.err.println(decodeCore + " <-decoded");
//            decodeCore = decodeCore(stringLength);
//        }
    }
}
