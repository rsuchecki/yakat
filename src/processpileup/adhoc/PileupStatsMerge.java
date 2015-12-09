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
package processpileup.adhoc;

import agrparser.ArgParser;
import agrparser.Opt;
import agrparser.OptSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PileupStatsMerge {

    private final String TOOL_NAME;
    private int READER_BUFFER_SIZE = 8192;
    private final int HELP_WIDTH = 200;

    public PileupStatsMerge(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
//        ArrayList<String> SAMPLE_NAMES = (ArrayList<String>) optSet.getOpt("n").getValues();
//        int minCoverageThreshold = (int) optSet.getOpt("c").getValueOrDefault();
//        int maxCoverageThreshold = (int) optSet.getOpt("C").getValueOrDefault();
        String PILEUP_FILE = (String) optSet.getOpt("p").getValueOrDefault();
//        Opt opt = optSet.getOpt("s");
//        String SNP_FILE = (String) opt.getValueOrDefault();
//        PAIRWISE_OVERLAP = !optSet.getOpt("P").isUsed();
        pileupStatsMerge(PILEUP_FILE);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('n', "sample-names", "list of sample names, order must correspond to mpileup input content").setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE).setRequired(true));
//        optSet.addOpt(new Opt('c', "min-coverage-per-sample", "Minimum coverage required for a sample to be processed/considered", 2, 1, 1000000));
//        optSet.addOpt(new Opt('C', "max-coverage-per-sample", "Maximum coverage allowed for a sample to be processed/considered", 1000, 1, 1000000));
        optSet.addOpt(new Opt('p', "pileup-file", "Input (m)pileup file, alternatively use stdin", 1));
//        optSet.addOpt(new Opt('s', "snp-file", "List of snps which will be used to interrogate offspring coverage at (parent) SNP positions", 1));
//        optSet.addOpt(new Opt('P', "no-pairwise-overlap", "Do not calculate pairwise overlap - this should speed things up"));

        return optSet;
    }

    private void pileupStatsMerge(String inputFile) {

        BufferedReader content = null;
        try {
            if (inputFile == null || inputFile.equals("-")) { //READSTDIN
                content = new BufferedReader(new InputStreamReader(System.in), READER_BUFFER_SIZE);
            } else if (inputFile.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile));
                content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                content = new BufferedReader(new FileReader(new File(inputFile)), READER_BUFFER_SIZE);
            }
            String line;
            long[][] matrix = new long[19][260];

            int index = 0;
            while ((line = content.readLine()) != null) {
                if(line.startsWith("cov")) {
                    index = 0;
                } else if(line.isEmpty()) {
                    
                } else {
                    String[] toks = line.split("\t");
                    for (int i = 1; i < toks.length; i++) {
                        matrix[i-1][index] += Long.parseLong(toks[i]);
//                        System.out.println("adding "+Long.parseLong(toks[i])+ " to ["+index+"]["+i+"]");
                    }
                    index++;
                }
            }
            
            long maxLen = String.valueOf(getMaxOfMatrix(matrix)).length()+1;
            for (int i = 0; i < matrix[0].length; i++) {
                System.out.printf("%" + 4 + "d", i+1);
                for (int j = 0; j < matrix.length; j++) {
                    System.out.printf("%" + maxLen + "d", matrix[j][i]);
                }
                System.out.println();
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found exception: " + ex.getMessage(), TOOL_NAME);
            System.exit(1);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        } finally {
            try {
                if (content != null) {
                    content.close();
                }
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
    }
    
    private long getMaxOfMatrix(long[][] matrix) {
        long max = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] > max) {
                    max = matrix[i][j];
                }
            }
        }
        return max;
    }

}
