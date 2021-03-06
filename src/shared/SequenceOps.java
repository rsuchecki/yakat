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

package shared;

/**
 * Generic methods for handling nucleotide strings.
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class SequenceOps {
    
    public static char complement(char c) {
        switch (c) {
            case 'A':
                c = 'T';
                break;
            case 'T':
                c = 'A';
                break;
            case 'C':
                c = 'G';
                break;
            case 'G':
                c = 'C';
                break;
            case 'N':
                c = 'N';
                break;
            case 'a':
                c = 't';
                break;
            case 't':
                c = 'a';
                break;
            case 'c':
                c = 'g';
                break;
            case 'g':
                c = 'c';
                break;
            case 'n':
                c = 'n';
                break;
            case  ' ':
                c = ' ';
            default: {
//                c = Character.toUpperCase(c);
//                System.err.println("Error reverse-complementing, fallen over char: " + c);
//                System.exit(1);
                break;
            }
        }
        return c;
    }
    
    public static String getReverseComplementString(CharSequence sequenceChars) {
        StringBuilder sb = new StringBuilder(sequenceChars.length());
        for (int i = sequenceChars.length() - 1; i > -1; --i) {
            sb.append(complement(sequenceChars.charAt(i)));
        }
        return sb.toString();
    }
    
    public static String getCanonical(String sequence) {
        String rc = getReverseComplementString(sequence);
        if(rc.compareTo(sequence) < 0) {
            return rc;
        } else {
            return sequence;
        }
    }
    
    public static CharSequence getReverseComplement(CharSequence sequenceChars) {
        StringBuilder sb = new StringBuilder(sequenceChars.length());
        for (int i = sequenceChars.length() - 1; i > -1; --i) {
            sb.append(complement(sequenceChars.charAt(i)));
        }
        return sb;
    }
    
    /**
     * 
     * @param sequence
     * @param from inclusive
     * @param to inclusive
     * @return 
     */
    public static boolean isCanonical(CharSequence sequence, int from, int to) {        
        int j = to;
        for (int i = from; i <= to ; i++) {
            int diff = sequence.charAt(i) - complement(sequence.charAt(j--));
            if(diff < 0) {
                return true;
            } else if (diff > 0) {
                return false;
            }            
        }
        return true; //is canonical when palindromic ()
    }
    
    public static boolean isCanonical(CharSequence sequence) {        
        return isCanonical(sequence, 0, sequence.length()-1);
    }
    
    public static CharSequence translate(CharSequence sequence) {        
        StringBuilder sb = new StringBuilder(sequence.length());
        CodonDictionary cd = new CodonDictionary();
        for (int i = 0; i < sequence.length()-2; i+=3) {
            sb.append(cd.translate(sequence.subSequence(i, i+3)));
        }
        return sb;
    }

    
}
