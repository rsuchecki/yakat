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

import java.util.ArrayList;

/**
 *
 * @author rad
 */
public class MsaSeqPair {

    private final MsaSequence s1;
    private final MsaSequence s2;
    private final ArrayList<Snp> snpsList;

    public MsaSeqPair(MsaSequence s1, MsaSequence s2) {
        if (s1.getId().compareTo(s2.getId()) < 0) {
            this.s1 = s1;
            this.s2 = s2;
        } else {
            this.s1 = s2;
            this.s2 = s1;
        }
        snpsList = new ArrayList<>();
    }

    public MsaSequence getS1() {
        return s1;
    }

    public MsaSequence getS2() {
        return s2;
    }

    public void addSnp(Snp snp) {
        char n1 = getS1().getSequenceString().charAt(snp.getPositionZeroIndexed());
        char n2 = getS2().getSequenceString().charAt(snp.getPositionZeroIndexed());
        if (n1 == '-' || n2 == '-') {
            if (!snpsList.isEmpty()) {
                Snp last = snpsList.get(snpsList.size() - 1);
                if (last instanceof Indel && snp.getPositionZeroIndexed() == last.getPositionZeroIndexed() + 1) {
                    last.incrementLength();
                } else {
                    snpsList.add(new Indel(s1, s2, snp.getPositionZeroIndexed()));
                }
            } else {
                snpsList.add(new Indel(s1, s2, snp.getPositionZeroIndexed()));
            }
        } else {
            snpsList.add(snp);
        }
//        this.snps++;
    }

//    public int getSnpsCount() {
//        return snpsList.size();
//    }
    
    /**
     * TODO - reduce gap extension penalty - currently weighted == gap open == mismatch
     * @return 
     */
    private int getSnpsWeight() {
        int weight = 0;
        for (Snp snp : snpsList) {
            int length = snp.getLength();
            weight += length > 0 ? length : 1;
        }
        return weight;
    }

    /**
     * Identity based on unpadded length of longer sequence
     *
     * @return
     */
    public double getMinIdentity() {
//        System.err.println(getSnpsWeight()+" "+s1.getUnpaddedLength()+" "+s2.getUnpaddedLength());
//        System.err.println(s1.getId() + " " + getIdentity1() + " " + getIdentity2() + " " + s2.getId());
        return Math.min(getIdentity1(), getIdentity2());
    }

    /**
     * Identity based on unpadded length of sequence 1
     *
     * @return
     */
    private double getIdentity1() {
        return 1 - (double) getSnpsWeight()/ s1.getUnpaddedLength();
    }

    /**
     * Identity based on unpadded length of sequence 2
     *
     * @return
     */
    private double getIdentity2() {
        return 1 - (double) getSnpsWeight()/ s2.getUnpaddedLength();
    }

}
