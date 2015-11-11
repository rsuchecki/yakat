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
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerStrings {

    private String kmer1String;
    private String kmer2String;

    public PairMerStrings(PairMer mer, int k) {
        if (!mer.isInvalid() && mer.getStoredCount() == 2) {

            String decodedCore = mer.decodeCore(k - 1);
            StringBuilder sb1 = new StringBuilder().append(mer.getClipLeft());
            sb1.append(decodedCore);
            this.kmer1String = sb1.toString();

            //Reconstruct k-mer 2
            StringBuilder sb2 = new StringBuilder(decodedCore);
            sb2.append(mer.getClipRight());
            this.kmer2String = sb2.toString();
        }
    }

    public String getKmer1String() {
        return kmer1String;
    }

    public String getKmer2String() {
        return kmer2String;
    }

    public String getOtherCoreOfKmer1() {
        return kmer1String.substring(0, kmer1String.length() - 1);
    }

    public String getOtherCoreOfKmer2() {
        return kmer2String.substring(1);
    }

}
