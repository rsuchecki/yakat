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
public class PairMer3LongEncoded extends PairMer implements Comparable<PairMer3LongEncoded> {

    private long kmerCoreBits1;
    private long kmerCoreBits2;
    private long kmerCoreBits3;

    /**
     * Proper constructor
     *
     * @param leftClip
     * @param core
     * @param rightClip
     */
    public PairMer3LongEncoded(char leftClip, String core, char rightClip, int freq) {
        addFirstKmer(leftClip, core, rightClip, freq);
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

    /**
     * Does not generate a complete PairMer, just the core, for Set/Maps lookups
     *
     * @param kmerCoreOnly
     */
    public PairMer3LongEncoded(String kmerCoreOnly) {
        encodeCore(SequenceOps.getCanonical(kmerCoreOnly));
    }

    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMer3LongEncoded) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        return CoreCoder.computeHash(getBitFields());
    }

    @Override
    public int compareTo(PairMer3LongEncoded anotherKmer) {
        return CoreCoder.compareCores(getBitFields(), anotherKmer.getBitFields());
    }

    public long[] getBitFields() {
        long bitsArray[] = {kmerCoreBits1, kmerCoreBits2, kmerCoreBits3};
        return bitsArray;
    }

    @Override
    public String decodeCore(int coreLength) {
        long bitsArray[] = {kmerCoreBits1, kmerCoreBits2, kmerCoreBits3};
        return CoreCoder.decodeCore(coreLength, bitsArray);
    }

    private void encodeCore(String kmerCoreOnly) {
        long[] encodeCoreLong = CoreCoder.encodeCoreLongArray(kmerCoreOnly);
        if (encodeCoreLong.length != 3) {
            Reporter.report("[BUG?]", " 3*long values expected from core encoding", getClass().getSimpleName());
        } else {
            kmerCoreBits1 = encodeCoreLong[0];
            kmerCoreBits2 = encodeCoreLong[1];
            kmerCoreBits3 = encodeCoreLong[2];
        }
//        String decodeCore = decodeCore(kmerCoreOnly.length());
//        if(!decodeCore.equals(kmerCoreOnly)) {
//            System.err.println("error");
//            System.err.println(kmerCoreOnly);
//            System.err.println(decodeCore);
//        }
    }

}
