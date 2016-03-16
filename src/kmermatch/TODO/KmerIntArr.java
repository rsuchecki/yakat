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
package kmermatch.TODO;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerIntArr implements Kmer, Comparable<KmerIntArr> {
    //Object overhead 8 B
    private int[] kmerBitsArray;  //12B + len*4 Bytes 

    public KmerIntArr(CharSequence sequence, int from, int to) {
        encodeKmerCanonical(sequence, from, to);
    }

    
//    @Override
    public String decodeKmer(int k) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int compareTo(KmerIntArr o) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public int[] getKmerBitsArray() {
        return kmerBitsArray;
    }
    
    /**
     * Encode kmer core in it's two possible representations (forward, rev-comp) but store canonical representation only
     * (decided by lex order)
     *
     * @param sequence
     * @param from inclusive
     * @param to inclusive
     * @return true if stored forward, false if RC
     */
    private boolean encodeKmerCanonical(CharSequence sequence, int from, int to) {
//        System.err.println("Encoding seq len=" + kmerString.length());

        int INT_LENGTH = 30; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
//        int stringLength = kmerString.length();
        int stringLength = to - from + 1;
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
//        int position = 0;
        int position = from;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        int[] kmerCoreBitsArray = new int[intsNeeded];
        int[] kmerCoreBitsArrayRC = new int[intsNeeded];
        int currentIntRc = intsNeeded - 1;
//        char[] kmerCharArray = kmerString.toCharArray();
        int positionInIntRc = 0;//positionInInt;

        while (position <= to) {
            while ((positionInInt < INT_LENGTH) && (position <= to)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                if (positionInIntRc >= INT_LENGTH) {
                    --currentIntRc;
                    positionInIntRc = 0;
                }
                switch (sequence.charAt(position)) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        //11 in Reverse-Complement                        
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        //10 in Reverse-Complement                        
                        positionInIntRc++;
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        //01 in Reverse-Complement                        
                        kmerCoreBitsArrayRC[currentIntRc] ^= 1 << positionInIntRc++;
                        positionInIntRc++;

                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        //00 in Reverse-Complement                        
                        positionInIntRc += 2;
                        break;
                    default:
                        System.err.println("Failed ecoding kmerstring to int array....");
                        System.err.println("Offending char: " + sequence.charAt(position));
                        System.err.println("in " + sequence);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }

        return storeCanonical(kmerCoreBitsArray, kmerCoreBitsArrayRC, stringLength);
    }
    
      /**
     *
     * @param forward
     * @param reverseComplement
     * @return true if stored in forward orientation, false if in reverse-complement
     */
    private boolean storeCanonical(int[] forward, int[] reverseComplement, int k) {

        for (int i = 0; i < forward.length; i++) {
            if (forward[i] < reverseComplement[i]) {
                break;
            } else if (forward[i] > reverseComplement[i]) {
                kmerBitsArray = reverseComplement;
                return false;
            }
        }
        kmerBitsArray = forward;
        return true;
    }

    @Override
    public String getKmerString() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
