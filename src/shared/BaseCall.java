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
package shared;

/**
 *
 * @author rad
 */
public class BaseCall {

    private final Character allele1;
    private final Character allele2;
    private final double freq1;
    private final double freq2;
    private final int cov1;
    private final int cov2;

    public BaseCall(Character allele1, Character allele2, double freq1, double freq2, int cov1, int cov2) {
        this.allele1 = allele1;
        this.allele2 = allele2;
        this.freq1 = freq1;
        this.freq2 = freq2;
        this.cov1 = cov1;
        this.cov2 = cov2;
    }

    public BaseCall() {
        this.allele1 = null;
        this.allele2 = null;
        this.freq1 = 0;
        this.freq2 = 0;
        this.cov1 = 0;
        this.cov2 = 0;
    }

    /**
     * TODO - not dealing with coverage and freq
     *
     * @param iupac
     */
    public BaseCall(char iupac) {
        switch (iupac) {
            case 'W':
                allele1 = 'A';
                allele2 = 'T';
                break;
            case 'S':
                allele1 = 'C';
                allele2 = 'G';
                break;
            case 'M':
                allele1 = 'A';
                allele2 = 'C';
                break;
            case 'K':
                allele1 = 'G';
                allele2 = 'T';
                break;
            case 'R':
                allele1 = 'A';
                allele2 = 'G';
                break;
            case 'Y':
                allele1 = 'C';
                allele2 = 'T';
                break;
            case 'A':
                allele1 = 'A';
                allele2 = null;
                break;
            case 'C':
                allele1 = 'C';
                allele2 = null;
                break;
            case 'G':
                allele1 = 'G';
                allele2 = null;
                break;
            case 'T':
                allele1 = 'T';
                allele2 = null;
                break;
            default:
                allele1 = null;
                allele2 = null;
                Reporter.report("[WARNING]", "Unable to convert IUPAC call " + iupac + " into a BaseCall object", this.getClass().getSimpleName());
        }
        this.freq1 = 0;
        this.freq2 = 0;
        this.cov1 = 0;
        this.cov2 = 0;
    }

    public boolean isCalled() {
        return allele1 !=null || allele2 != null;
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

    public String getCoverageString() {
        return getCov1() + "/" + getCov2();
    }

    public String getFreqString() {
        return getFreq1() + "/" + getFreq2();
    }

    public CharSequence getEvidence() {
        String EMPTY = "\u2205"; 
        String TIMES = "\u00D7";
        StringBuilder sb = new StringBuilder();
        if (getCov1() == 0) {
            sb.append(EMPTY);
        } else {
            sb.append(getCov1()).append(TIMES).append(getFreq1());
        }
        sb.append("/");
        if (getCov2() == 0) {
            sb.append(EMPTY);
        } else {
            sb.append(getCov2()).append(TIMES).append(getFreq2());
        }                
        return sb;
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

    public double getFreq1() {
        return freq1;
    }

    public String getFreq1String() {
        String f = Double.toString(getFreq1());
        return f.replaceAll("\\.0+$", "");
    }
    

    public double getFreq2() {
        return freq2;
    }
    
    public String getFreq2String() {
        String f = Double.toString(getFreq2());
        return f.replaceAll("\\.0+$", "");
    }

    public int getCov1() {
        return cov1;
    }

    public int getCov2() {
        return cov2;
    }

}
