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
package freqmers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import shared.LabelledInputBuffer;
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class CallerConsumer implements Runnable {

    private final BlockingQueue<LabelledInputBuffer> inputQueue;
    private final String TOOL_NAME;
    private final ArrayList<KmerFilter> snpFilters;
//    private final ConcurrentSkipListMap<CharSequence, KmerLink> map;
    private final HashMap<String, ArrayList<KmerLink>> map;
//    private final OptSet optSet;
    private final ArrayList<String> samples;

//    ConcurrentHashMap<String, PerSampleBuffer> sampleToBufferMap;
//    ConcurrentHashMap<String, BlockingQueue<PerSampleBuffer>> sampleToQueueMap;
    public CallerConsumer(BlockingQueue<LabelledInputBuffer> inputQueue, String TOOL_NAME, ArrayList<String> samples,
            HashMap<String, ArrayList<KmerLink>> map, ArrayList<KmerFilter> snpFilters) {
        this.inputQueue = inputQueue;
        this.samples = samples;
        this.map = map;
        this.snpFilters = snpFilters;
//        this.optSet = optSet;
        this.TOOL_NAME = TOOL_NAME;
//        ONLY_COUNT = optSet.getOpt("only-count").getOptFlag();
//        MIN_LENGTH_READ = (int) optSet.getOpt("r").getValueOrDefault();
//        sampleToBufferMap = new ConcurrentHashMap<>(keyMap.getSamplesTotal() * 2);
//        sampleToQueueMap = keyMap.getSampleToQueueMap();
    }

    @Override
    public void run() {

        try {
            LabelledInputBuffer laeblledBuffer = null;
            String previous = null;
            while (!(laeblledBuffer = inputQueue.take()).getData().isEmpty()) {
                if (!samples.contains(laeblledBuffer.getLabel())) {  //FIRST OR NEW SAMPLE (not seen before)
                    if (!samples.isEmpty()) { //NEW SAMPLE
                        Reporter.report("[INFO]", "Collecting k-mer frequencies for " + previous + " and reseting mer-counters", TOOL_NAME);
                        for (KmerFilter snpFilter : snpFilters) {
                            snpFilter.collectStatsAndResetMers(samples.get(samples.size() - 1), TOOL_NAME);
                        }
                    }
//                    Reporter.report("[INFO]", "Current sample: " + list.getLabel(), TOOL_NAME);
                    previous = laeblledBuffer.getLabel();
                    samples.add(laeblledBuffer.getLabel());
                }
                ArrayList<String[]> data = laeblledBuffer.getData();
                for (String[] toks : data) {
                    ArrayList<KmerLink> kmerLinks = map.get(SequenceOps.getCanonical(toks[0]));
                    if (kmerLinks != null) {
//                        if(kmerLinks.size()>1) {
//                            int f=0;
//                        }
                        for (KmerLink kmerLink : kmerLinks) {
                            if (kmerLink != null) {
                                boolean setMer = kmerLink.setMer((short) Math.min(Integer.parseInt(toks[1]), Short.MAX_VALUE)); 
                                if (!setMer) {
                                    Reporter.report("[ERROR]", "Unable to set k-mer link to sequence, possible reason: duplicate k-mer in an input set, k-mer: " + toks[0], TOOL_NAME);
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
            Reporter.report("[INFO]", "Collecting k-mer frequencies for " + previous + " and reseting mer-counters", TOOL_NAME);
            for (KmerFilter kmerFilter : snpFilters) {
                kmerFilter.collectStatsAndResetMers(samples.get(samples.size() - 1), TOOL_NAME);
            }
            Reporter.report("[INFO]", "Finished assigning k-mers to sequences", TOOL_NAME);

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
