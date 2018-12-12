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
package shared;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public final class CodonDictionary {

    private final String[] CODON_DICTIONARY = {"Ala/A	GCT", "Ala/A	GCC", "Ala/A	GCA", "Ala/A	GCG", "Leu/L	TTA", "Leu/L	TTG", "Leu/L	CTT", "Leu/L	CTC", "Leu/L	CTA", "Leu/L	CTG", "Arg/R	CGT", "Arg/R	CGC", "Arg/R	CGA", "Arg/R	CGG", "Arg/R	AGA", "Arg/R	AGG", "Lys/K	AAA", "Lys/K	AAG", "Asn/N	AAT", "Asn/N	AAC", "Met/M	ATG", "Asp/D	GAT", "Asp/D	GAC", "Phe/F	TTT", "Phe/F	TTC", "Cys/C	TGT", "Cys/C	TGC", "Pro/P	CCT", "Pro/P	CCC", "Pro/P	CCA", "Pro/P	CCG", "Gln/Q	CAA", "Gln/Q	CAG", "Ser/S	TCT", "Ser/S	TCC", "Ser/S	TCA", "Ser/S	TCG", "Ser/S	AGT", "Ser/S	AGC", "Glu/E	GAA", "Glu/E	GAG", "Thr/T	ACT", "Thr/T	ACC", "Thr/T	ACA", "Thr/T	ACG", "Gly/G	GGT", "Gly/G	GGC", "Gly/G	GGA", "Gly/G	GGG", "Trp/W	TGG", "His/H	CAT", "His/H	CAC", "Tyr/Y	TAT", "Tyr/Y	TAC", "Ile/I	ATT", "Ile/I	ATC", "Ile/I	ATA", "Val/V	GTT", "Val/V	GTC", "Val/V	GTA", "Val/V	GTG", "STOP/*	TAA", "STOP/*	TGA", "STOP/*	TAG"};
    private final HashMap<CharSequence, CharSequence> codonDictionary;

    public CodonDictionary(String codonsDictionaryFileName) {
        codonDictionary = codonDictionaryReader(codonsDictionaryFileName);
    }

    public CodonDictionary() {
        codonDictionary = codonDictionaryReader(CODON_DICTIONARY);

    }

    public CharSequence translate(CharSequence codon) {
        if (codon.length() == 3) {
            CharSequence aa = codonDictionary.get(((String) codon).toUpperCase());     
            if(codon == "TAA") {
                System.err.println(codon+" -> "+codonDictionary.get(((String) codon).toUpperCase()));
            }
            return aa == null ? "X" : aa;
        } else {
            return "?";
        }
    }

    private HashMap<CharSequence, CharSequence> codonDictionaryReader(String[] codonsDictionary) {
        HashMap<CharSequence, CharSequence> dictionary = new HashMap<>();
        for (String inputLine : codonsDictionary) {
            String splits[] = inputLine.split(" |\t");
            String key = splits[1].toUpperCase();
            dictionary.put(key, splits[0].replaceFirst(".*/", ""));
        }
        return dictionary;
    }

    private HashMap<CharSequence, CharSequence> codonDictionaryReader(String codonsDictionaryFileName) {
        HashMap<CharSequence, CharSequence> dictionary = new HashMap<>();
        File newFile = new File(codonsDictionaryFileName);
        BufferedReader myData = null;
        try {
            String inputLine;
            myData = new BufferedReader(new FileReader(newFile));
            while ((inputLine = myData.readLine()) != null) {
                String splits[] = inputLine.split(" |\t");
                String key = splits[1].toUpperCase();
                CharSequence previous = dictionary.put(key, splits[0].replaceFirst(".*/", ""));

                if (previous != null) {
                    System.err.println("Fatal error: non unique entry: " + key + " in " + codonsDictionaryFileName);
                    System.exit(5);
                }
            }
        } catch (FileNotFoundException ex) {
            System.err.println("File not found exception!\n\t" + newFile.getName());
//            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (myData != null) {
                    myData.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        return dictionary;
    }
}
