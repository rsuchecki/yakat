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
package vsearchprocess;

/**
 *
 * @author rad
 */
public class MsaSeqPair {
    private MsaSequence s1;
    private MsaSequence s2;
    private int snps = 0;

    public MsaSeqPair(MsaSequence s1, MsaSequence s2) {
        if (s1.getId().compareTo(s2.getId()) < 0) {
            this.s1 = s1;
            this.s2 = s2;
        } else {
            this.s1 = s2;
            this.s2 = s1;
        }
    }

    public MsaSequence getS1() {
        return s1;
    }

    public MsaSequence getS2() {
        return s2;
    }

    public void addSnp() {
        this.snps++;
    }

    public int getSnps() {
        return snps;
    }
    
    /**
     * Identity based on unpadded length of longer sequence 
     * 
     * @return 
     */
    public double getMinIdentity() {
        return Math.min(getIdentity1(), getIdentity2());    
    }
    
    /**
     * Identity based on unpadded length of sequence 1
     * @return 
     */
    public double getIdentity1() {
        return 1 - (double)snps/s1.getUnpaddedLength();
    }
    
    /**
     * Identity based on unpadded length of sequence 2
     * @return 
     */
    public double getIdentity2() {
        return 1 - (double)snps/s2.getUnpaddedLength();
    }

}
