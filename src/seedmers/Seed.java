/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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
package seedmers;

import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Seed {

    private final String id;
    private final CharSequence sequence;
    private final int snpPosition;
    private short mersA[];
    private short mersC[];
    private short mersG[];
    private short mersT[];

    public Seed(String id, CharSequence sequence, int k) {
        this.id = id;
        this.sequence = sequence;
        this.snpPosition = k-1;
    }

    public String getId() {
        return id;
    }

    public CharSequence getSequence() {
        return sequence;
    }

    public int getSnpPosition() {
        return snpPosition;
    }

    public boolean setMer(int position, short value, String kmer) {
        char allele = kmer.charAt(getSnpPosition() - position + 1);
        short mers[];
        switch (allele) {
            case 'A':
                mers = mersA;
                break;
            case 'C':
                mers = mersC;
                break;
            case 'G':
                mers = mersG;
                break;
            case 'T':
                mers = mersT;
                break;
            default: {
                Reporter.report("[WARNING]", "Unrecognized based in kmer "+kmer+ " - failed setting mer count...", this.getClass().getSimpleName());
                return false;
            }
        }        
        if (mers[position] == 0) {
            mers[position] = value;
            return true;
        }
        return false;
    }
}
