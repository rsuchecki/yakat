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


/**
 * Proof of concept structure to hold up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerIntArrEncodedWithRef extends PairMerIntArrEncoded {

    private final PairMer next;

    public PairMerIntArrEncodedWithRef(CharSequence sequence, int from, int to, boolean frontClip, int freq, PairMer next) throws NonACGTException {
        super(sequence, from, to, frontClip, freq);
        this.next = next;
    }

    public PairMer getNext() {
        return next;
    }
    
    @Override
    public PairMer getNextPairMer() {
        return next;
    }

    
    
}
