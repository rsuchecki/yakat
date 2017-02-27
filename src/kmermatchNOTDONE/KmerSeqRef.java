/*
 * Copyright 2016 rad.
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

import shared.SequenceOps;

/**
 *
 * @author rad
 */
public class KmerSeqRef implements Kmer, Comparable<Kmer> {
    CharSequence kmerSequence;

    
    public KmerSeqRef(CharSequence sequence, int from, int to, boolean canonical) {
        if(canonical) {            
            kmerSequence = SequenceOps.getCanonical(sequence.subSequence(from, to).toString());            
        } else {
            kmerSequence = sequence.subSequence(from, to);            
        }
    }

//    KmerSeqRef(CharSequence kmerSequence, int from, int to) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }

    @Override
    public int compareTo(Kmer another) {
        return getKmerString().compareTo(another.getKmerString());
    }
    
    @Override
    public String getKmerString() {
        return kmerSequence instanceof String ? (String)kmerSequence : kmerSequence.toString();
    }
    
}
