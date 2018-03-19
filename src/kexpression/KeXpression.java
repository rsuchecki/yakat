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
package kexpression;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import shared.Reporter;
import shared.StdRedirect;
import snpmers.SnpMers;

/**
 * Mostly a placeholder,
 *
 * @author rad
 */
public class KeXpression {

    private final static String DELIMITER = "\t";
    private final static double READ_LEN = 100;
    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int READER_BUFFER_SIZE = 8192;

    public KeXpression(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
//        new StdRedirect(optSet, TOOL_NAME);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

//        if (optSet.getOpt("F").isUsed()) {
//        }
//        int minTotal = (int) optSet.getOpt("min-k-mer-frequency-sum").getValueOrDefault();
//        int minMinor = (int) optSet.getOpt("min-k-mer-frequency-minor").getValueOrDefault();
//        double minCoverage = (double) optSet.getOpt("min-snp-coverage").getValueOrDefault();
//        double maxError = (double) optSet.getOpt("max-coverage-error").getValueOrDefault();
        GenesMap genesMap = parseNotQuiteGff((String) optSet.getOpt("q").getValueOrDefault());
        computeFeaturesLengths(genesMap);
        computeTPMs(genesMap, (String) optSet.getOpt("c").getValueOrDefault());
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("Currently: Calculate TPMs, start by processing transcript+exon info to establish feature lengths");

        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('q', "not-quite-gff", "Gene lengths will be computed from GFF-derived transcript and exon coordinates", 1));
        optSet.addOpt(new Opt('c', "read-counts-matrix", "These will be converted to TPMs", 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    private GenesMap parseNotQuiteGff(String fileName) {
        GenesMap genesMap = new GenesMap();

        BufferedReader bufferdReader = null;
        try {
            String inputLine;
            if (fileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(fileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(fileName)), READER_BUFFER_SIZE);
            }
            while ((inputLine = bufferdReader.readLine()) != null) {
                String line = inputLine.trim();
                String[] toks = line.split(DELIMITER);
                if (toks[0].equals("transcript")) {
                    genesMap.addTranscript(toks[3], new Transcript(toks[4], Integer.parseInt(toks[1]), Integer.parseInt(toks[2])));
                } else if (toks[0].equals("exon")) {
                    genesMap.addExon(toks[3], toks[4], new Exon(Integer.parseInt(toks[1]), Integer.parseInt(toks[2])));
                }
//                System.out.println(line);
            }
            Reporter.report("[INFO]", genesMap.getGenes().size() + " genes parsed.", TOOL_NAME);
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + fileName, TOOL_NAME);
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bufferdReader != null) {
                    bufferdReader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        return genesMap;
    }

    private void computeFeaturesLengths(GenesMap genesMap) {
        HashMap<String, Gene> genes = genesMap.getGenes();
        for (Map.Entry<String, Gene> entry : genes.entrySet()) {
            String geneId = entry.getKey();
            Gene g = entry.getValue();
            HashMap<String, Transcript> transcripts = g.getTranscripts();
            //RECORD OUTER COORDINATES OF A FEATURE (GENE)
            int from = Integer.MAX_VALUE;
            int to = 0;
            for (Map.Entry<String, Transcript> tEntry : transcripts.entrySet()) {
//                String transcriptId = tEntry.getKey();
                Transcript t = tEntry.getValue();
                from = t.getFrom() < from ? t.getFrom() : from;
                to = t.getTo() > to ? t.getTo() : to;
            }
            //CHECK WHICH OF THE BASES WITHIN AND INCLUDING THE OUTER COORDINATES ARE EXONIC
            int length = 0;
            for (int i = from; i <= to; i++) {
                scanExons:
                for (Map.Entry<String, Transcript> tEntry : transcripts.entrySet()) {
                    Transcript t = tEntry.getValue();
                    ArrayList<Exon> exons = t.getExons();
                    for (Exon exon : exons) {
                        if (i >= exon.getFrom() && i <= exon.getTo()) {
                            length++;
                            break scanExons;
                        }
                    }
                }
            }
            g.setLength(length);
//            System.out.println(geneId + DELIMITER + length);
        }
        Reporter.report("[INFO]", "Finished computing feature lengths...", TOOL_NAME);
    }

    private void computeTPMs(GenesMap genesMap, String fileName) {
        BufferedReader bufferdReader = null;
        try {
            String inputLine;
            if (fileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(fileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(fileName)), READER_BUFFER_SIZE);
            }
            String header[] = null;
            //COLLECT COUNTS            
            while ((inputLine = bufferdReader.readLine()) != null) {
                String line = inputLine.trim();
                String[] toks = line.split(DELIMITER);
                if (toks[0].equals("GeneId")) { //PROCESS HEADER
                    header = toks;
                } else {
                    for (int i = 1; i < toks.length; i++) {

                        Gene gene = genesMap.getGenes().get(toks[0]);
                        if (gene == null) {
                            System.err.println("Failed to find gene! " + toks[0]);
                        } else {
                            gene.addMappedReads(header[i], Integer.parseInt(toks[i]));
                        }
                    }
                }
//                System.out.println(line);
            }
            TreeSet geneIds = new TreeSet(genesMap.getGenes().keySet());
            //COMPUTE Ts
            double perSampleTs[] = new double[header.length];
            for (int i = 1; i < header.length; i++) {
                String sampleId = header[i];
                Iterator<String> it = geneIds.iterator();
                while (it.hasNext()) {
                    String geneId = it.next();
                    perSampleTs[i] += genesMap.getMappedReadsPerGene(geneId, sampleId) * READ_LEN / genesMap.getGeneLength(geneId);
                }
            }
            Reporter.report("[INFO]", "Finished computing Ts for normalising...", TOOL_NAME);
            //            
            StringBuilder line = new StringBuilder(header[0]);
            for (int i = 1; i < header.length; i++) {
                line.append(DELIMITER).append(header[i]).append("[").append(perSampleTs[i]).append("]");
            }
            System.out.println(line);
            //COMPUTE TPMs
            double mil = 1000000;
            Iterator<String> it = geneIds.iterator();
            while (it.hasNext()) {
                String geneId = it.next();
                line = new StringBuilder(geneId);
                for (int i = 1; i < header.length; i++) {
                    String sampleId = header[i];
//                    System.out.println(genesMap.getMappedReadsPerGene(geneId, sampleId)+"*100*"+mil
//                    +"/"+genesMap.getGeneLength(geneId)+"*"+perSampleTs[i]);
                    double tpm = genesMap.getMappedReadsPerGene(geneId, sampleId) * READ_LEN * mil / genesMap.getGeneLength(geneId) * perSampleTs[i];
                    line.append(DELIMITER).append(tpm);
                }
                System.out.println(line);
            }

//            Reporter.report("[INFO]", genesMap.getGenes().size() + " genes parsed.", TOOL_NAME);
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + fileName, TOOL_NAME);
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bufferdReader != null) {
                    bufferdReader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    private class GenesMap {

        private HashMap<String, Gene> genes;

        public GenesMap() {
            genes = new HashMap<>();
        }

        public void addTranscript(String geneId, Transcript t) {
//            String geneId = t.getId().replaceFirst("\\.[0-9]+", "");
            Gene gene = genes.get(geneId);
            if (gene == null) {
                gene = new Gene(geneId);
                genes.put(geneId, gene);
            }
            gene.addTranscript(t);
        }

        public void addExon(String geneId, String transcriptId, Exon e) {
            genes.get(geneId).getTranscript(transcriptId).addExon(e);
        }

        public HashMap<String, Gene> getGenes() {
            return genes;
        }

        public int getGeneLength(String geneId) {
            return genes.get(geneId).getLength();
        }

        public int getMappedReadsPerGene(String geneId, String sampleId) {
            return genes.get(geneId).getMappedReads(sampleId);
        }
    }

    private class Gene {

        private HashMap<String, Double> tpmsMap;
        private HashMap<String, Integer> mappedReadMap;
        private HashMap<String, Transcript> transcripts;
        private final String id;
        private int length;

        public Gene(String id) {
            this.id = id;
            transcripts = new HashMap<>();
            mappedReadMap = new HashMap<>();
        }

        public void addTranscript(Transcript t) {
            transcripts.put(t.getId(), t);
        }

        public Transcript getTranscript(String transcriptId) {
            return transcripts.get(transcriptId);
        }

        public HashMap<String, Transcript> getTranscripts() {
            return transcripts;
        }

        public int getLength() {
            return length;
        }

        public void setLength(int length) {
            this.length = length;
        }

        public void addMappedReads(String sample, int value) {
            mappedReadMap.put(sample, value);
        }

//        public void addTpm(String sample, double value) {
//            tpmsMap.put(sample, value);
//        }
        public Integer getMappedReads(String sample) {
            Integer value = mappedReadMap.get(sample);
            if (value == null) {
                return 0;
            }
            return value;
        }

    }

    private class Transcript {

        private ArrayList<Exon> exons;
        private final int from;
        private final int to;
        private final String id;

        public Transcript(String id, int from, int to) {
            this.id = id;
            exons = new ArrayList<>();
            this.from = from;
            this.to = to;
        }

        public void addExon(Exon e) {
            exons.add(e);
        }

        public int getFrom() {
            return from;
        }

        public int getTo() {
            return to;
        }

        public String getId() {
            return id;
        }

        public ArrayList<Exon> getExons() {
            return exons;
        }

    }

    private class Exon {

        private final int from;
        private final int to;

        public Exon(int from, int to) {
            this.from = from;
            this.to = to;
        }

        public int getFrom() {
            return from;
        }

        public int getTo() {
            return to;
        }
    }
}
