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
package main;

import allorfs.AllOrfs;
import argparser.ParserTest;
import fastqmatchid.FastqMatchId;
import freqmers.FreqMers;
import gbssplit.SplitGBS;
import hmmerdoms.HmmerDomains;
import kexpression.KeXpression;
import kmerextender.CoreCoder;
import kmerextender.KmerExtender;
import kmermatch.KmerMatch;
//import kmermatchNOTDONE.KmerMatch;
import kmerger.KmerSetMerge;
import processpileup.ProcessPileup;
import processpileup.PileupStats;
import pileupstatsmerge.PileupStatsMerge;
import seedmers.SeedMers;
import snpmers.SnpMers;
import vsearchprocess.VsearchClustersCaller;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class MerUtensils {

    public static void main(String[] args) {
        new MerUtensils(args);
    }

    public MerUtensils(String[] args) {
        String name = this.getClass().getSimpleName().toLowerCase();
        if (args.length != 0 && args[0].matches("(k)?extend(er)?")) {
            new KmerExtender(args, name, "kextend");
//        } else if (args.length != 0 && args[0].matches("test(k)?extend(er)?")) {
//            new CoreCoder(args, name, "testkextender");
        } else if (args.length != 0 && args[0].matches("(k)?match(er)?")) {
            args[0] = "";
            new KmerMatch(args, name, "kmatch");
//            new KmerMatch(args, name, "kmatch");
        } else if (args.length != 0 && args[0].matches("(k)?merge(r)?")) {
            new KmerSetMerge(args, name, "kmerge");
        } else if (args.length != 0 && args[0].matches("split(er)?")) {
            new SplitGBS(args, name, "split");
        } else if (args.length != 0 && args[0].matches("idmatch(er)?")) {
            new FastqMatchId(args, name, "idmatch");
        } else if (args.length != 0 && args[0].matches("pileupstats")) {
            new PileupStats(args, name, "pileup");
        } else if (args.length != 0 && args[0].matches("(p)?mpileup")) {
            new ProcessPileup(args, name, "pmpileup");
        } else if (args.length != 0 && args[0].matches("pileupstatsmerge")) {
            new PileupStatsMerge(args, name, "pileupstatsmerge");
        } else if (args.length != 0 && args[0].matches("test")) {
            args[0] = "";
            new ParserTest(this.getClass().getSimpleName(), "test", args);
        } else if (args.length != 0 && args[0].matches("sn(i)?pmers")) {
            new SnpMers(args, name, "snpmers");
        } else if (args.length != 0 && args[0].matches("s(eed)?mers")) {
            new SeedMers(args, name, "seedmers");
        } else if (args.length != 0 && args[0].matches("f(req)?mers")) {
            new FreqMers(args, name, "freqmers");
        } else if (args.length != 0 && args[0].matches("h(a)?mmerdom(ain)?s")) {
            new HmmerDomains(args, name, "hmmerdoms");
        } else if (args.length != 0 && args[0].matches("(all)?orfs")) {
            new AllOrfs(args, name, "allorfs");
        } else if (args.length != 0 && args[0].matches("vclust(ers)?")) {
//            new (args, name, "vclustersvars");
            new VsearchClustersCaller(args, name, "vclust");
        } else if (args.length != 0 && args[0].matches("(k)?express(ion)?")) {
            new KeXpression(args, name, "kexpress");
        } else if (args.length != 0 && args[0].matches("(v|ver)(sion)?")) {
            Package aPackage = this.getClass().getPackage();
            String version = aPackage.getImplementationVersion();
//            String build = aPackage.getImplementationTitle();
//            System.out.println(version + " " + build);
            System.out.println(version);
        } else {
            printHelp();
        }
    }

    private void printHelp() {
        System.out.println();
        String version = this.getClass().getPackage().getImplementationVersion();
        System.out.println("java -jar " + this.getClass().getSimpleName().toLowerCase() + "-" + version + ".jar <command> ");
        System.out.println("Commands/modules:");
        System.out.println("  k-mer based utils");
        System.out.println("    kextend       : extend k-mers to unambiguous contigs (optionally extend input \"seed\" sequences only)");
        System.out.println("    snpmers       : given parental SNPs and the corresponding FASTA sequences (e.g. from NIKS), call offspring genotypes by overlapping their k-mers with parental SNP sequences ");
        System.out.println("    freqmers      : given a set of sequences and set(s) of k-mers, report k-mer coverage and frequency for the input sequences");
        System.out.println("    kmatch        : match/filter/bait FAST(A|Q) sequences based on contained k-mers (or lack thereof)");        
        System.out.println("    seedmers      : [PROTOTYPE] given seed seequences interrogare sets of k-mers to genotype presumed mutations at positions k bases from the seed edges");
        System.out.println("    kmerge        : [DEPRECATED - use KMC tools for k-mer set operations] given sorted input, merge k-mer sets summing frequencies if available");
//        System.out.println("   kmatch     : a.k.a bait");
        System.out.println();
        System.out.println("  FASTQ processing utils:");
        System.out.println("    idmatch       : match FASTQ records by id");
        System.out.println("    split         : split FASTQ GBS/ddRAD reads by barcodes, trim barcodes and adapters");
        System.out.println();
        System.out.println("  MPILEUP processing utils:");
        System.out.println("    pileupstats   : extract some stats from (m)pileup");
        System.out.println("    pmpileup      : count and call bases from (m)pileup");
        System.out.println();
        System.out.println("  HMMER, VSEARCH output post-processing utils:");
        System.out.println("    hmmerdoms     : Process and group domain-hits from HMMER");
        System.out.println("    vclusters     : call variants from VSEARCH clustering msa output");
        System.out.println();
        System.out.println("  Miscellaneous utils:");
        System.out.println("    allorfs       : Identify and extract all (longest) ORFs from a genome");
        System.out.println("    kexpress      : [UNTESTED] Calculate expresion values (TPM) from read counts and a (not quite) GFF description of features ");
//        System.out.println("   mmers      : count (and analyse?) k-mers in multiple input sets ");
        System.out.println();
        System.out.println("  Info:");
        System.out.println("    version       : print the version and exit");

//        String s = "Currently k-mer frequency is not taken into consideration, so use of a dedicated k-mer counting program, "
//                + "such as KMC or Jellyfish is recommended. It is best to exclude low frequency k-mers before passing "
//                + "the list of k-mers to KmerExtender. For smaller jobs FASTA or FASTQ input may suffice.";
//        System.out.println(Reporter.wrapString(s, 145));
        System.out.println();
    }
}
