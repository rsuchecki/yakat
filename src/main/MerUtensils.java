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

import agrparser.ParserTest;
import kmerextender.KmerExtender;
import kmermatch.KmerMatch;
import kmersetmerge.KmerSetMerge;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class MerUtensils {

    public static void main(String[] args) {
        new MerUtensils(args);
    }

    public MerUtensils(String[] args) {
        if (args.length != 0 && args[0].matches("(k)?extend(er)?")) {
            args[0] = "";
            new KmerExtender(args);
        } else if (args.length != 0 && args[0].matches("(k)?match(er)?")) {
            args[0] = "";
            new KmerMatch(args);
        } else if (args.length != 0 && args[0].matches("(k)?merge(r)?")) {       
            new KmerSetMerge(args, this.getClass().getSimpleName());
        } else if (args.length != 0 && args[0].matches("test")) {
            args[0] = "";
            new ParserTest(this.getClass().getSimpleName(), "test", args);
        } else {
            printHelp();
        }
    }

    private void printHelp() {
        System.out.println();
        System.out.println("java -jar " + this.getClass().getSimpleName() + ".jar <command> ");
        System.out.println("Commands:");
        System.out.println("   extend  : ");
        System.out.println("   match   : ");
        System.out.println("   merge   : given sorted input, merge k-mer sets summing frequencies if available ");
//        String s = "Currently k-mer frequency is not taken into consideration, so use of a dedicated k-mer counting program, "
//                + "such as KMC or Jellyfish is recommended. It is best to exclude low frequency k-mers before passing "
//                + "the list of k-mers to KmerExtender. For smaller jobs FASTA or FASTQ input may suffice.";
//        System.out.println(Reporter.wrapString(s, 145));
        System.out.println();
    }
}
