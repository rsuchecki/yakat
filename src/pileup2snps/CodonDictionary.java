/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pileup2snps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public final class CodonDictionary {

    private final String[] CODON_DICTIONARY = {"Ala/A	GCT", "Ala/A	GCC", "Ala/A	GCA", "Ala/A	GCG", "Leu/L	TTA", "Leu/L	TTG", "Leu/L	CTT", "Leu/L	CTC", "Leu/L	CTA", "Leu/L	CTG", "Arg/R	CGT", "Arg/R	CGC", "Arg/R	CGA", "Arg/R	CGG", "Arg/R	AGA", "Arg/R	AGG", "Lys/K	AAA", "Lys/K	AAG", "Asn/N	AAT", "Asn/N	AAC", "Met/M	ATG", "Asp/D	GAT", "Asp/D	GAC", "Phe/F	TTT", "Phe/F	TTC", "Cys/C	TGT", "Cys/C	TGC", "Pro/P	CCT", "Pro/P	CCC", "Pro/P	CCA", "Pro/P	CCG", "Gln/Q	CAA", "Gln/Q	CAG", "Ser/S	TCT", "Ser/S	TCC", "Ser/S	TCA", "Ser/S	TCG", "Ser/S	AGT", "Ser/S	AGC", "Glu/E	GAA", "Glu/E	GAG", "Thr/T	ACT", "Thr/T	ACC", "Thr/T	ACA", "Thr/T	ACG", "Gly/G	GGT", "Gly/G	GGC", "Gly/G	GGA", "Gly/G	GGG", "Trp/W	TGG", "His/H	CAT", "His/H	CAC", "Tyr/Y	TAT", "Tyr/Y	TAC", "Ile/I	ATT", "Ile/I	ATC", "Ile/I	ATA", "Val/V	GTT", "Val/V	GTC", "Val/V	GTA", "Val/V	GTG", "STOP/*	TAA", "STOP/*	TGA", "STOP/*	TAG"};
    private final HashMap<String, String> codonDictionary;

    public CodonDictionary(String codonsDictionaryFileName) {
        codonDictionary = codonDictionaryReader(codonsDictionaryFileName);
    }

    public CodonDictionary() {
        codonDictionary = codonDictionaryReader(CODON_DICTIONARY);

    }

    public String translate(String codon) {
        if(codon.length() == 3)
        return codonDictionary.get(codon.toUpperCase());
        else
            return "?";
    }

    public HashMap<String, String> codonDictionaryReader(String[] codonsDictionary) {
        HashMap<String, String> dictionary = new HashMap<>();
        for (String inputLine : codonsDictionary) {
            String splits[] = inputLine.split(" |\t");
            String key = splits[1].toUpperCase();
            dictionary.put(key, splits[0].replaceFirst(".*/", ""));
        }
        return dictionary;
    }

    private HashMap<String, String> codonDictionaryReader(String codonsDictionaryFileName) {
        HashMap<String, String> dictionary = new HashMap<>();
        File newFile = new File(codonsDictionaryFileName);
        BufferedReader myData = null;
        try {
            String inputLine;
            myData = new BufferedReader(new FileReader(newFile));
            while ((inputLine = myData.readLine()) != null) {
                String splits[] = inputLine.split(" |\t");
                String key = splits[1].toUpperCase();
                String previous = dictionary.put(key, splits[0].replaceFirst(".*/", ""));

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
