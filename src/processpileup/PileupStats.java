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

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
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
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class PileupStats {

    private final String TOOL_NAME;
    boolean IGNORE_SPLICING = true;
    private boolean SAMPLE_VS_REF = false;
    private boolean PAIRWISE_OVERLAP = true;
//    private String INPUT_FILE;
    private int READER_BUFFER_SIZE = 8192;
//    private ArrayList<String> SAMPLE_NAMES;
    private final int HELP_WIDTH = 200;

    public PileupStats(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        ArrayList<String> SAMPLE_NAMES = (ArrayList<String>) optSet.getOpt("n").getValues();
        int minCoverageThreshold = (int) optSet.getOpt("c").getValueOrDefault();
        int maxCoverageThreshold = (int) optSet.getOpt("C").getValueOrDefault();
        String PILEUP_FILE = (String) optSet.getOpt("p").getValueOrDefault();
        Opt opt = optSet.getOpt("s");
        String SNP_FILE = (String) opt.getValueOrDefault();
        PAIRWISE_OVERLAP = !optSet.getOpt("P").isUsed();
        pileupStats(PILEUP_FILE, SAMPLE_NAMES, minCoverageThreshold, maxCoverageThreshold);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('n', "sample-names", "list of sample names, order must correspond to mpileup input content").setMinValueArgs(2).setMaxValueArgs(Integer.MAX_VALUE).setRequired(true));
        optSet.addOpt(new Opt('c', "min-coverage-per-sample", "Minimum coverage required for a sample to be processed/considered", 2, 1, 1000000));
        optSet.addOpt(new Opt('C', "max-coverage-per-sample", "Maximum coverage allowed for a sample to be processed/considered", 1000, 1, 1000000));
        optSet.addOpt(new Opt('p', "pileup-file", "Input (m)pileup file, alternatively use stdin", 1));
        optSet.addOpt(new Opt('s', "snp-file", "List of snps which will be used to interrogate offspring coverage at (parent) SNP positions", 1));
        optSet.addOpt(new Opt('P', "no-pairwise-overlap", "Do not calculate pairwise overlap - this should speed things up"));

        return optSet;
    }

    private void pileupStats(String inputFile, ArrayList<String> sampleIds, int minCoverageThreshold, int maxCoverageThreshold) {

        int[][] pairwiseOverlapAtCoverage;

        int[] cumulativeCoverage = new int[maxCoverageThreshold +2 ];
        if (PAIRWISE_OVERLAP) {
            pairwiseOverlapAtCoverage = new int[sampleIds.size()][sampleIds.size()];
        } else {
            pairwiseOverlapAtCoverage = null;
        }
        int[][] samplesAtMinCoverage = new int[21][sampleIds.size() + 1];
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
            while ((line = content.readLine()) != null && !line.isEmpty()) {
                String[] toks = line.split("\t");
                int[] record = new int[sampleIds.size()];
                int r = 0;
                int cumCoverage = 0;
                //tos 0,1,2 are ref, position, refbase                                
                for (int i = 3; i < toks.length; i += 3) {
                    try {
                        int coverage = Integer.parseInt(toks[i]);
                        record[r++] = coverage;
                        cumCoverage += coverage;
                    } catch (ArrayIndexOutOfBoundsException e) {
                        Reporter.report("[FATAL]", "Array index out of bounds - likely cause: mismatch between samples given and pileup file content", TOOL_NAME);
                        System.exit(1);
                    }
                }
                if (cumCoverage > maxCoverageThreshold+1) {
                    cumulativeCoverage[maxCoverageThreshold]++;
                } else {
                    cumulativeCoverage[cumCoverage]++;
                }
                if (PAIRWISE_OVERLAP) {
                    for (int i = 0; i < pairwiseOverlapAtCoverage.length; i++) {
                        for (int j = 0; j < pairwiseOverlapAtCoverage[0].length; j++) {
                            if (record[i] >= minCoverageThreshold && record[j] >= minCoverageThreshold && record[i] <= maxCoverageThreshold && record[j] <= maxCoverageThreshold) {
                                pairwiseOverlapAtCoverage[j][i]++;
                            }
                        }
                    }
                }
                //samples @ coverage distribution
                for (int c = 2; c < samplesAtMinCoverage.length; c++) {
                    int recordsAtOrOverCoverage = 0;
                    for (int i = 0; i < record.length; i++) {
                        if (record[i] >= c) {
                            recordsAtOrOverCoverage++;
                        }
                    }
                    samplesAtMinCoverage[c][recordsAtOrOverCoverage]++;
                }
            }
            if (PAIRWISE_OVERLAP) {
                int maxLen = Math.max(String.valueOf(getMaxOfMatrix(pairwiseOverlapAtCoverage)).length() + 1, getMaxNameLen(sampleIds)) + 1;
                System.out.println();
                System.out.printf("%" + maxLen + "s", "");
                for (String id : sampleIds) {
                    System.out.printf("%" + maxLen + "s", id);
                }
                System.out.println();
                for (int i = 0; i < pairwiseOverlapAtCoverage.length; i++) {
                    System.out.printf("%" + maxLen + "s", sampleIds.get(i));
                    for (int j = 0; j < pairwiseOverlapAtCoverage[0].length; j++) {
                        System.out.printf("%" + maxLen + "d", pairwiseOverlapAtCoverage[i][j]);
                    }
                    System.out.println();
                }
            }
            int maxLen = Math.max(String.valueOf(getMaxOfMatrix(samplesAtMinCoverage)).length() + 1, ("cov>=" + samplesAtMinCoverage.length).length() + 1);
            System.out.println();
            System.out.printf("%" + String.valueOf(sampleIds.size()).length() + "s", "");
            for (int i = 0; i < samplesAtMinCoverage.length; i++) {
                if (i < 2) {
//                    System.out.printf("%" + maxLen + "s", "");
                } else {
                    System.out.printf("%" + maxLen + "s", "cov>=" + i);
                }
            }
            System.out.println();
            for (int i = 0; i < samplesAtMinCoverage[0].length; i++) {
                System.out.printf("%" + String.valueOf(sampleIds.size()).length() + "d", i);
                for (int j = 2; j < samplesAtMinCoverage.length; j++) {
                    System.out.printf("%" + maxLen + "d", samplesAtMinCoverage[j][i]);
                }
                System.out.println();
            }
//            System.out.println();
//            for (int i = 1; i < cumulativeCoverage.length; i++) {
//                if (i <= maxCoverageThreshold) {
//                    System.err.println(i + "\t" + cumulativeCoverage[i]);
//                } else {
//                    System.err.println("more\t" + cumulativeCoverage[i]);                    
//                }
//            }

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
        //    private HashMap<String,PileupPosition> parseSnpsTable(String fName) {
    //        
    //    }

    private int getMaxOfMatrix(int[][] matrix) {
        int max = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] > max) {
                    max = matrix[i][j];
                }
            }
        }
        return max;
    }

    private int getMaxNameLen(ArrayList<String> names) {
        int max = 0;
        for (String name : names) {
            if (name.length() > max) {
                max = name.length();
            }
        }
        return max;
    }

}
