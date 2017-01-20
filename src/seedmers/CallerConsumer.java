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
 *//*
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
package seedmers;

import argparser.OptSet;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.regex.Pattern;
import shared.LabelledInputBuffer;
import shared.Message;
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class CallerConsumer implements Runnable {

    private final BlockingQueue<LabelledInputBuffer> inputQueue;
    private final String TOOL_NAME;
    private final HashMap<String, ArrayList<AltSeedLink>> map;
//    private final OptSet optSet;
    private final ArrayList<String> samples;
//    private final int minTotal;
//    private final int minMinor;
//    private final double minCoverage;
//    private final double maxError;
    ArrayList<Seed> seeds;

//    ConcurrentHashMap<String, PerSampleBuffer> sampleToBufferMap;
//    ConcurrentHashMap<String, BlockingQueue<PerSampleBuffer>> sampleToQueueMap;
    public CallerConsumer(BlockingQueue<LabelledInputBuffer> inputQueue, String TOOL_NAME, ArrayList<String> samples,
            HashMap<String, ArrayList<AltSeedLink>> map, ArrayList<Seed> seeds) {
//            HashMap<String, ArrayList<AltSeedLink>> map, int minTotal, int minMinor, double minKmers, double maxError) {
        this.inputQueue = inputQueue;
        this.samples = samples;
        this.map = map;
//        this.snpFilters = altSeedLinks;
//        this.optSet = optSet;
        this.TOOL_NAME = TOOL_NAME;
        this.seeds = seeds;
//        ONLY_COUNT = optSet.getOpt("only-count").getOptFlag();
//        MIN_LENGTH_READ = (int) optSet.getOpt("r").getValueOrDefault();
//        sampleToBufferMap = new ConcurrentHashMap<>(keyMap.getSamplesTotal() * 2);
//        sampleToQueueMap = keyMap.getSampleToQueueMap();
//        this.minTotal = minTotal;
//        this.minMinor = minMinor;
//        this.minCoverage = minKmers;
//        this.maxError = maxError;
    }

    @Override
    public void run() {

//         char labels[] = {'A','C','G','T'};
//        StringBuilder header = new StringBuilder("ID");
//        for (char label : labels) {
//            header.append("\t");
//            header.append("#kmers_").append(label);
//            header.append("\t");
//            header.append("median_cov_").append(label);
//            header.append("\t");
//            header.append("cov_ratio_").append(label);
//            
//        }
//        System.out.println(header.toString());
        System.out.println("ID\t#kmers\tmedian_cov\tcov_ratio\tCall");
        
        try {
            LabelledInputBuffer laeblledBuffer = null;
            String previous = null;
            while (!(laeblledBuffer = inputQueue.take()).getData().isEmpty()) {
                if (!samples.contains(laeblledBuffer.getLabel())) {  //FIRST OR NEW SAMPLE (not seen before)
                    if (!samples.isEmpty()) { //NEW SAMPLE
                        Reporter.report("[INFO]", "Calling bases for " + previous + " and reseting mer-counters", TOOL_NAME);
                        for (Seed seed : seeds) {
                            seed.callBaseAndResetMers(samples.get(samples.size() - 1), TOOL_NAME);
                        }
                    }
//                    Reporter.report("[INFO]", "Current sample: " + list.getLabel(), TOOL_NAME);
                    previous = laeblledBuffer.getLabel();
                    samples.add(laeblledBuffer.getLabel());
                }
                ArrayList<String[]> data = laeblledBuffer.getData();
                for (String[] toks : data) {
                    String canonical = SequenceOps.getCanonical(toks[0]);
//                    System.err.println("Current "+toks[0]+" canonical: "+canonical);
                    ArrayList<AltSeedLink> altSeedLinks = map.get(canonical);
                    if (altSeedLinks != null) {
//                        if(kmerLinks.size()>1) {
//                            int f=0;
//                        }
                        for (AltSeedLink altSeedLink : altSeedLinks) {
                            if (altSeedLink != null) {
//                                System.err.println("Setting "+toks[0]);
                                boolean setMer = altSeedLink.getParent().setMer(altSeedLink, (short) Math.min(Integer.parseInt(toks[1]), Short.MAX_VALUE));
                                if (!setMer) {
                                    Reporter.report("[ERROR]", "Unable to set k-mer link to SNP, possible reason: duplicate k-mer in an input set, k-mer: " + toks[0], TOOL_NAME);
                                }
//                    SnpFilter snpFilter = kmerLink.getSnpFilter();
//                    System.err.println(kmerLink.getParentSequence().getId()+"\t"+snpFilter.getSnpPosition0()+"\t"+toks[1]);
                            }
                        }
                    }
                }

            }
//            inputQueue.put(new LabelledInputBuffer(null, new ArrayList())); //inform other threads
            //PROCESS LAST SAMPLE RESULTS
            Reporter.report("[INFO]", "Calling bases for " + previous + " and reseting mer-counters", TOOL_NAME);
            for (Seed seed : seeds) {
                seed.callBaseAndResetMers(samples.get(samples.size() - 1), TOOL_NAME);
            }
            Reporter.report("[INFO]", "Finished assigning k-mer frequencies to candidate variants", TOOL_NAME);

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
