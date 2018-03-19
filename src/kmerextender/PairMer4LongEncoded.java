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
 * Storing PairMer core in 3 long fields
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMer4LongEncoded extends PairMer implements Comparable<PairMer4LongEncoded> {

    private long kmerCoreBits1;
    private long kmerCoreBits2;
    private long kmerCoreBits3;
    private long kmerCoreBits4;

    /**
     * Proper constructor
     *
     * @param sequence
     * @param from
     * @param to
     * @param frontClip
     * @param freq
     */
    public PairMer4LongEncoded(CharSequence sequence, int from, int to, boolean frontClip, int freq) {
        addFirstKmer(sequence, from, to, frontClip, freq);
    }

    public final void addFirstKmer(CharSequence sequence, int from, int to, boolean frontClip, int freq) {
        if (getStoredCount() == 0) {        //If this is the first of the two k-mers that could be stored
            int coreStart = frontClip ? from + 1 : from;
            int coreEnd = frontClip ? to : to - 1;
            if (SequenceOps.isCanonical(sequence.subSequence(coreStart, coreEnd+1))) {
                encodeCore(sequence.subSequence(coreStart, coreEnd+1));
                if (frontClip) {
                    setClipLeft(sequence.charAt(from));
                } else {
                    setClipRight(sequence.charAt(to));
                }
            } else {
                encodeCore(SequenceOps.getReverseComplement(sequence.subSequence(coreStart, coreEnd+1)));
                if (frontClip) {
                    setClipRight(SequenceOps.complement(sequence.charAt(from)));
                } else {
                    setClipLeft(SequenceOps.complement(sequence.charAt(to)));
                }
            }
            incrementStoredCount(hasLeftClip(), freq);
        } else {
            Reporter.report("[BUG?]", "Only the first k-mer in a PairMer can be added using addFirstKmer()!!!", getClass().getSimpleName());
        }
    }

    /**
     * Does not generate a complete PairMer, just the core, for Set/Map lookups
     *
     * @param kmerCoreOnly
     */
    public PairMer4LongEncoded(CharSequence kmerCoreOnly) {
        if (SequenceOps.isCanonical(kmerCoreOnly)) {
            encodeCore(kmerCoreOnly);
        } else {
            encodeCore(SequenceOps.getReverseComplement(kmerCoreOnly));
        }
    }

    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMer4LongEncoded) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        return CoreCoder.computeHash(getBitFields());
    }

    @Override
    public int compareTo(PairMer4LongEncoded anotherKmer) {
        return CoreCoder.compareCores(getBitFields(), anotherKmer.getBitFields());
    }

    public long[] getBitFields() {
        long bitsArray[] = {kmerCoreBits1, kmerCoreBits2, kmerCoreBits3, kmerCoreBits4};
        return bitsArray;
    }

    @Override
    public String decodeCore(int coreLength) {
        long bitsArray[] = {kmerCoreBits1, kmerCoreBits2, kmerCoreBits3, kmerCoreBits4};
        return CoreCoder.decodeCore(coreLength, bitsArray);
    }

    private void encodeCore(CharSequence kmerCoreOnly) {
        long[] encodeCoreLong = CoreCoder.encodeCoreLongArray(kmerCoreOnly);
        if (encodeCoreLong.length != 4) {
            Reporter.report("[BUG?]", " 4*long values expected from core encoding", getClass().getSimpleName());
        } else {
            kmerCoreBits1 = encodeCoreLong[0];
            kmerCoreBits2 = encodeCoreLong[1];
            kmerCoreBits3 = encodeCoreLong[2];
            kmerCoreBits4 = encodeCoreLong[3];
        }
//        //Sanity check
//        String decodeCore = decodeCore(kmerCoreOnly.length());
//        if(!decodeCore.equals(kmerCoreOnly.toString())) {
//            System.err.println("error");
//            System.err.println(kmerCoreOnly);
//            System.err.println(decodeCore);
//        }
    }

}
