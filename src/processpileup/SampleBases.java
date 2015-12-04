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
package processpileup;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SampleBases {

    private int aCount;
    private int cCount;
    private int gCount;
    private int tCount;
//    private short indelSize;
    private int indelCount;
    private int otherCount;
    private int pileupDepth;

    public SampleBases(String pileupRecord, String refBase, int pileupDepth, boolean IGNORE_SPLICING) {
        this.pileupDepth = pileupDepth;
        processPileupRecord(pileupRecord, refBase, IGNORE_SPLICING);
    }
    
    
    
    
    public SampleBases processPileupRecord(String pileupRecord, String refBase, boolean IGNORE_SPLICING) {
        if (refBase == null) {
            pileupRecord = pileupRecord.toUpperCase().replaceAll("\\*", "I").replaceAll("\\^.|\\$", "");;
        } else {
            pileupRecord = pileupRecord.toUpperCase().replaceAll("\\.|,", refBase).replaceAll("\\*", "I").replaceAll("\\^.|\\$", "");;
        }

        for (int i = 0; i < pileupRecord.length(); i++) {
            char base = pileupRecord.charAt(i);

            if ((base == '+' || base == '-') && i + 2 < pileupRecord.length()) {
                int numIndel = Integer.parseInt(pileupRecord.substring(i + 1, i + 2));
//                key = pileupRecord.substring(i, i + numIndel + 2);
                i += numIndel + 2;
                base = 'i';
            }
            if (IGNORE_SPLICING && (base == '<' || base == '>')) {
                continue;
            }
            switch (base) {
                case 'A':
                    aCount++;
                    break;
                case 'C':
                    cCount++;
                    break;
                case 'G':
                    gCount++;
                    break;
                case 'T':
                    tCount++;
                    break;
                case 'i':
                    indelCount++;
//                    indelSize =
                    break;
                default:
                    otherCount++;
                    break;
            }
        }
        return this;
    }
    
    

}
