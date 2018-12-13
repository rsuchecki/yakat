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

import shared.SequenceOps;
import shared.Reporter;

/**
 * Storing PairMer core in 2 long fields
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class PairMer1LongEncoded extends PairMer implements Comparable<PairMer1LongEncoded> {

    private long kmerCoreBits1;

    /**
     * Proper constructor
     *
     * @param sequence
     * @param from
     * @param to
     * @param frontClip
     * @param freq
     */
    public PairMer1LongEncoded(CharSequence sequence, int from, int to, boolean frontClip, int freq) {
        addFirstKmer(sequence, from, to, frontClip, freq);
    }

    

    /**
     * Does not generate a complete PairMer, just the core, for Set/Map lookups
     *
     * @param kmerCoreOnly
     */
    public PairMer1LongEncoded(CharSequence kmerCoreOnly) {
        if (SequenceOps.isCanonical(kmerCoreOnly)) {
            encodeCore(kmerCoreOnly);
        } else {
            encodeCore(SequenceOps.getReverseComplement(kmerCoreOnly));
        }
    }

    @Override
    public boolean equals(Object anotherKmer) {
        return compareTo((PairMer1LongEncoded) anotherKmer) == 0;
    }

    @Override
    public int hashCode() {
        return CoreCoder.computeHash(getBitFields());
    }

    @Override
    public int compareTo(PairMer1LongEncoded anotherKmer) {
        return CoreCoder.compareCores(getBitFields(), anotherKmer.getBitFields());
    }

    public long[] getBitFields() {
        long bitsArray[] = {kmerCoreBits1};
        return bitsArray;
    }

    @Override
    public String decodeCore(int coreLength) {
        long bitsArray[] = {kmerCoreBits1};
        return CoreCoder.decodeCore(coreLength, bitsArray);
    }

    @Override
    protected void encodeCore(CharSequence kmerCoreOnly) {
        long[] encodeCoreLong = CoreCoder.encodeCoreLongArray(kmerCoreOnly);
        if (encodeCoreLong.length != 1) {
            Reporter.report("[BUG?]", " 1*long value expected from core encoding", getClass().getSimpleName());
        } else {
            kmerCoreBits1 = encodeCoreLong[0];
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
