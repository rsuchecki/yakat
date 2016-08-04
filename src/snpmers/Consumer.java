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
package snpmers;

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
public class Consumer implements Runnable {

    private final BlockingQueue<LabelledInputBuffer> inputQueue;
    private final String TOOL_NAME;
    private final ArrayList<Message> finalMessages;
    private final ArrayList<SnpFilter> snpFilters;
//    private final ConcurrentSkipListMap<CharSequence, KmerLink> map;
    private final HashMap<CharSequence, KmerLink> map;
    private final OptSet optSet;
    private final ArrayList<String> samples;
    
//    ConcurrentHashMap<String, PerSampleBuffer> sampleToBufferMap;
//    ConcurrentHashMap<String, BlockingQueue<PerSampleBuffer>> sampleToQueueMap;
    public Consumer(BlockingQueue<LabelledInputBuffer> inputQueue, String TOOL_NAME, ArrayList<String> samples,
        HashMap<CharSequence, KmerLink> map, ArrayList<SnpFilter> snpFilters,
        OptSet optSet, ArrayList<Message> finalMessages) {
        this.inputQueue = inputQueue;
        this.samples = samples;
        this.map = map;
        this.snpFilters = snpFilters;
        this.optSet = optSet;
        this.TOOL_NAME = TOOL_NAME;
//        ONLY_COUNT = optSet.getOpt("only-count").getOptFlag();
//        MIN_LENGTH_READ = (int) optSet.getOpt("r").getValueOrDefault();
        this.finalMessages = finalMessages;
//        sampleToBufferMap = new ConcurrentHashMap<>(keyMap.getSamplesTotal() * 2);
//        sampleToQueueMap = keyMap.getSampleToQueueMap();
    }

    @Override
    public void run() {
        int minTotal = (int) optSet.getOpt("min-k-mer-frequency-sum").getValueOrDefault();
        int minMinor = (int) optSet.getOpt("min-k-mer-frequency-minor").getValueOrDefault();
        int minKmers = (int) optSet.getOpt("min-overlapping-k-mers").getValueOrDefault();
//        ArrayList<String> samples = new ArrayList<>();
        try {
            LabelledInputBuffer list = null;
            
            while (!(list = inputQueue.take()).getData().isEmpty()) {
                if (!samples.contains(list.getLabel())) {  //FIRST OR NEW SAMPLE
                    if (!samples.isEmpty()) { //NEW SAMPLE
                        Reporter.report("[INFO]", "Calling bases and reseting mer-counters", TOOL_NAME);
                        for (SnpFilter snpFilter : snpFilters) {
                            snpFilter.callBaseAndResetMers(samples.get(samples.size() - 1), minTotal, minMinor, minKmers, TOOL_NAME);
                        }
                    }
                    Reporter.report("[INFO]", "Current sample: " + list.getLabel(), TOOL_NAME);
                    samples.add(list.getLabel());
                }
                ArrayList<String[]> data = list.getData();
                for (String[] toks : data) {
                    KmerLink kmerLink = map.get(toks[0]);
                    if (kmerLink != null) {
                        boolean setMer = kmerLink.setMer(Short.parseShort(toks[1]));
                        if (!setMer) {
                            System.err.println("mer not set");
                        }
//                    SnpFilter snpFilter = kmerLink.getSnpFilter();
//                    System.err.println(kmerLink.getParentSequence().getId()+"\t"+snpFilter.getSnpPosition()+"\t"+toks[1]);
                    }
                }

            }
//            inputQueue.put(new LabelledInputBuffer(null, new ArrayList())); //inform other threads
            //PROCESS LAST SAMPLE RESULTS
            Reporter.report("[INFO]", "Calling bases and reseting mer-counters", TOOL_NAME);
            for (SnpFilter snpFilter : snpFilters) {
                snpFilter.callBaseAndResetMers(samples.get(samples.size() - 1), minTotal, minMinor, minKmers, TOOL_NAME);
            }
            Reporter.report("[INFO]", "Finished assigning k-mer frequencies to SNPs", TOOL_NAME);
//            if (lines > 0) {
//                String name = Thread.currentThread().getName();
//                String message = "[" + name + "] " + NumberFormat.getNumberInstance().format(lines)
//                    + " records processed, no matching barcode in " + NumberFormat.getNumberInstance().format(noBarcodeMatch);
//                if (OUTPUT_PSTI_STARTING_ONLY) {
//                    message += ", non-PstI start detected in " + NumberFormat.getNumberInstance().format(notPstIStart);
//                }
//                finalMessages.add(new Message(Message.Level.INFO, message, TOOL_NAME));
//                if (TRIM_ADAPTERS) {
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Trimmed total " + NumberFormat.getNumberInstance().format(trimmedMspIcount + trimmedBothCutSitesInPair + trimmedPstIcount) + " fragments", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Read-through detection with barcode and adapter trimming: ", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified MspI+3' Adapter in R1 and PstI+barcode in R2 in " + NumberFormat.getNumberInstance().format(trimmedBothCutSitesInPair) + " pairs", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified PstI+barcode in " + NumberFormat.getNumberInstance().format(trimmedPstIcount) + " R2 reads", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified MspI+3' Adapter in " + NumberFormat.getNumberInstance().format(trimmedMspIcount) + " R1 reads", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Length filtering breakdown: ", TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(pairUnderLenSum) + " pairs under combined length " + MIN_LENGTH_PAIR_SUM, TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(pairedReadUnderLen) + " paired reads under length " + MIN_LENGTH_PAIR_EACH, TOOL_NAME));
//                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(singleUnderLength) + " single reads under length " + MIN_LENGTH_READ, TOOL_NAME));
//                }
//            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
