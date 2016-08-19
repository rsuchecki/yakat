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
package snpmers;

import shared.Reporter;

/**
 *
 * @author rad
 */
public class BaseCall {

    private Character allele1;
    private Character allele2;

    /**
     *
     * @param allele1
     * @param allele2
     */
    public BaseCall(Character allele1, Character allele2) {
        this.allele1 = allele1;
        this.allele2 = allele2;
    }

    public String getCallString() {
        if (allele1 == null && allele2 == null) {
            return "N";
        }
        if (allele1 == null) {
            return allele2.toString();
        }
        if (allele2 == null) {
            return allele1.toString();
        }
        return allele1 + "/" + allele2;
    }

    public String getCallAB(char parent1call, char parent2call) {
        if (allele1 == null && allele2 == null) {
            return "NN";
        }
        if (allele1 == null) {
            if (allele2 == parent1call) {
                return "AA";
            }
            if (allele2 == parent2call) {
                return "BB";
            }
            System.err.println("BUG1 in AB call " + this.getClass());
            return "nn";
        }
        if (allele2 == null) {
            if (allele1 == parent1call) {
                return "AA";
            }
            if (allele1 == parent2call) {
                return "BB";
            }
            System.err.println("BUG2 in AB call " + this.getClass());
            return "nn";
        }
        if ((allele1 == parent1call && allele2 == parent2call) || (allele1 == parent2call && allele2 == parent1call)) {
            return "AB";
        }
        System.err.println("BUG3 in AB call " + this.getClass());
        return "nn";
    }

    public char getCallIUPAC() {
        if (allele1 == null && allele2 == null) {
            return 'N';
        }
        if (allele1 == null) {
            return allele2;
        }
        if (allele2 == null) {
            return allele1;
        }
        if (allele1 == '-' || allele2 == '-') {
            return 'n';
        }
        return getIupacCode(allele1, allele2);
    }

    private char getIupacCode(char c1, char c2) {
        if (c1 > c2) {
            char tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        if (c1 == 'A' && c2 == 'T') {
            return 'W';
        } else if (c1 == 'C' && c2 == 'G') {
            return 'S';
        } else if (c1 == 'A' && c2 == 'C') {
            return 'M';
        } else if (c1 == 'G' && c2 == 'T') {
            return 'K';
        } else if (c1 == 'A' && c2 == 'G') {
            return 'R';
        } else if (c1 == 'C' && c2 == 'T') {
            return 'Y';
        } else {
            Reporter.report("[FATAL] ", "Error calling IUPAC code for: " + c1 + "," + c2, this.getClass().getSimpleName());
            return '0';
        }
    }
}
