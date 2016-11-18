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
package kmermatchNOTDONE;

import kmerextender.*;
import java.util.Iterator;
import java.util.concurrent.ConcurrentSkipListMap;
import shared.SequenceOps;

/**
 * Wrapper around a concurrent collection (ConcurrentSkipListMap, but other options possible) of PairMer objects. After
 * the Map is populated, use .purge()
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerMap extends shared.MerMap {

    private final ConcurrentSkipListMap<Kmer, Kmer> kmersSkipListMap;

    private final int MAX_1LONG_ENCODE = 32; //2bits per nucl, signed long so should be 31, but can use sign bit if lex ordering not needed, so 32 allowed 
    private final int MAX_2LONG_ENCODE = 64;
    private final int MAX_3LONG_ENCODE = 96;
    private final int MAX_4LONG_ENCODE = 128;
    private final int MAX_5LONG_ENCODE = 160;

    /**
     * Instantiate the Map
     */
    public KmerMap() {
        kmersSkipListMap = new ConcurrentSkipListMap<>();
    }

    public void addToKmersMap(CharSequence sequence, int from, int to, boolean cannonical) {
        Kmer kmer;
//        if (kmerString.length() - 1 <= MAX_1LONG_ENCODE) {
//            kmer = new PairMer1LongEncoded(splitMer);
//        } else if (kmerString.length() - 1 <= MAX_2LONG_ENCODE) {
//            kmer = new PairMer2LongEncoded(splitMer);
//        } else if (kmerString.length() - 1 <= MAX_3LONG_ENCODE) {
//            kmer = new PairMer3LongEncoded(splitMer);
//        } else if (kmerString.length() - 1 <= MAX_4LONG_ENCODE) {
//            kmer = new PairMer4LongEncoded(splitMer);
//        } else if (kmerString.length() - 1 <= MAX_5LONG_ENCODE) {
//            kmer = new PairMer5LongEncoded(splitMer);
//        } else {
//            kmer = new KmerIntArrIntCounter(kmerString, count);
//        }
        kmer = new KmerSeqRef(sequence, from, to, cannonical);

        //Atomic operation START
        Kmer previousStoredKmer = kmersSkipListMap.putIfAbsent(kmer, kmer);
        //Atomic operation END
        if (previousStoredKmer != null) {
//            previousStoredKmer.incrementStoredCount(count);
        }
    }

    public Kmer get(CharSequence kmerSequence, int from, int to, boolean canonical) {

//        if (k - 1 <= MAX_1LONG_ENCODE) {
//            return kmersSkipListMap.get(new PairMer1LongEncoded(core));
//        } else if (k - 1 <= MAX_2LONG_ENCODE) {
//            return kmersSkipListMap.get(new PairMer2LongEncoded(core));
//        } else if (k - 1 <= MAX_3LONG_ENCODE) {
//            return kmersSkipListMap.get(new PairMer3LongEncoded(core));
//        } else if (k - 1 <= MAX_4LONG_ENCODE) {
//            return kmersSkipListMap.get(new PairMer4LongEncoded(core));
//        } else if (k - 1 <= MAX_5LONG_ENCODE) {
//            return kmersSkipListMap.get(new PairMer5LongEncoded(core));
//        } else {
//            return kmersSkipListMap.get(new KmerIntArrIntCounter(core));
//        }
        return kmersSkipListMap.get(new KmerSeqRef(kmerSequence, from, to, canonical));
    }

    /**
     *
     * @return the underlying data structure
     */
    public ConcurrentSkipListMap getPairMersSkipListMap() {
        return kmersSkipListMap;
    }

}
