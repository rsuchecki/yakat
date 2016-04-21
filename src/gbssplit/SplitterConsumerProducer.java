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
package gbssplit;

import argparser.OptSet;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;
import shared.Message;
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SplitterConsumerProducer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
//    private final Integer MIN_KMER_FREQUENCY;
//    private final Integer MAX_KMER_FREQUENCY;
//    private final String OUT_LABEL;
    private final int OUT_BUFFER_SIZE;
    private final KeyMap keyMap;
    private final String TOOL_NAME;
    private boolean TRIM_BARCODE = true;
    private boolean TRIM_ADAPTERS = true;
    private boolean PSTI_STARTING_ONLY = true;
    private boolean ONLY_COUNT = false; //no output produced if true

//    private AtomicInteger PRODUCER_THREADS;
    private final int MIN_LENGTH_READ;
    private final int MIN_LENGTH_PAIR_SUM;
    private final int MIN_LENGTH_PAIR_EACH;
    private final String ADAPTER;
    private ArrayList<Message> finalMessages;

    public SplitterConsumerProducer(BlockingQueue<ArrayList<String>> inputQueue,
        KeyMap keyMap, String toolName, int OUT_BUFFER_SIZE, OptSet optSet,
        ArrayList<Message> finalMessages) {

        this.inputQueue = inputQueue;
        this.keyMap = keyMap;
        this.TOOL_NAME = toolName;
        this.OUT_BUFFER_SIZE = OUT_BUFFER_SIZE;

        TRIM_BARCODE = !optSet.getOpt("keep-barcodes").getOptFlag();
        TRIM_ADAPTERS = !optSet.getOpt("keep-adapters").getOptFlag();
        PSTI_STARTING_ONLY = !optSet.getOpt("keep-non-PstI-starting").getOptFlag();
        ONLY_COUNT = optSet.getOpt("only-count").getOptFlag();

        ADAPTER = (String) optSet.getOpt("A").getValueOrDefault();
        MIN_LENGTH_READ = (int) optSet.getOpt("r").getValueOrDefault();
        MIN_LENGTH_PAIR_EACH = (int) optSet.getOpt("e").getValueOrDefault();
        MIN_LENGTH_PAIR_SUM = (int) optSet.getOpt("s").getValueOrDefault();
        this.finalMessages = finalMessages;
    }

    @Override
    public void run() {
        try {
            ConcurrentHashMap<String, SampleBuffer> sampleToBufferMap = new ConcurrentHashMap<>(keyMap.getSamplesTotal() * 2);
            ConcurrentHashMap<String, BlockingQueue<SampleBuffer>> sampleToQueueMap = keyMap.getSampleToQueueMap();
            ArrayList<String> list;
            Pattern spliPattern = Pattern.compile("\t");
            long noBarcodeMatch = 0L;
            long notPstIStart = 0L;
            long trimmedMspIcount = 0L;
            long trimmedPstIcount = 0L;
            long trimmedBothCutSitesInPair = 0L;
            long pairUnderLenSum = 0L;
            long pairedReadUnderLen = 0L;
            long singleUnderLength = 0L;
            long lines = 0L;
            while (!(list = inputQueue.take()).isEmpty()) {

//            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    lines++;
//                    StringTokenizer tokenizer = new StringTokenizer(line);
                    String toks[] = spliPattern.split(line);
//                    String id = tokenizer.nextToken();
//                    ConcurrentHashMap<String, String> barcodes = keyMap.getBarcodes(id.substring(1).replaceAll(":.*", ""));
                    ConcurrentHashMap<String, String> barcodes = keyMap.getBarcodes(toks[0].substring(1).replaceAll(":.*", ""));

                    if (barcodes != null) {
                        int matchingBarcodes = 0;
                        for (Map.Entry<String, String> entrySet : barcodes.entrySet()) {
                            String barcode = entrySet.getKey();
                            if (toks[1].startsWith(barcode)) {
                                matchingBarcodes++;
                                if (PSTI_STARTING_ONLY && !toks[1].startsWith(barcode + "TGCAG")) {
                                    notPstIStart++;
                                    break;
                                }
                                String sample = entrySet.getValue();
//                                ArrayList<String> bufferList = sampleToBufferMap.get(sample);
                                SampleBuffer bufferList = sampleToBufferMap.get(sample);
                                if (bufferList == null) { //nothing yet for this sample
                                    bufferList = new SampleBuffer(sample, null, OUT_BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);

                                }
                                if (bufferList.size() == OUT_BUFFER_SIZE) { //buffer full, place on queue to write and start a new buffer
                                    BlockingQueue<SampleBuffer> outputQueue = sampleToQueueMap.get(sample);
                                    keyMap.addToSampleCount(sample, OUT_BUFFER_SIZE);
                                    if (!ONLY_COUNT) {
                                        outputQueue.put(bufferList);
                                    }
                                    bufferList = new SampleBuffer(sample, null, OUT_BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);
                                }
                                //TRIM AND ASSESS RECORDS
                                StringBuilder builderR1 = new StringBuilder();
                                StringBuilder builderR2 = new StringBuilder();
                                builderR1.append(toks[0]); //id

                                if (TRIM_BARCODE) {
                                    toks[1] = toks[1].substring(barcode.length());
                                    toks[3] = toks[3].substring(barcode.length());
//                                    if (PSTI_STARTING_ONLY && !toks[1].startsWith("TGCAG")) {
//                                        continue;
//                                    }
//                                } else if (PSTI_STARTING_ONLY && !toks[1].startsWith(barcode + "TGCAG")) {

                                }
                                boolean trimmedMspI = false;
                                if (TRIM_ADAPTERS) {
                                    int trimFrom = toks[1].indexOf("CCG" + ADAPTER);
                                    if (trimFrom >= 0) {
//                                        System.err.println(toks[1]+" <--BEFORE l="+toks[1].length());
                                        toks[1] = toks[1].substring(0, trimFrom + 3);
//                                        System.err.println(toks[1]+" <--AFTER  l="+toks[1].length());
//                                        System.err.println(toks[3]+" <--BEFORE l="+toks[3].length());
                                        toks[3] = toks[3].substring(0, trimFrom + 3);
//                                        System.err.println(toks[3]+" <--AFTER  l="+toks[3].length());
                                        trimmedMspI = true;
                                    }
                                }
                                builderR1.append("\t").append(toks[1]); //sequence (possibly trimmed)
                                builderR1.append("\t").append("+"); //redundant id or '+
                                builderR1.append("\t").append(toks[3]); //qualline (possibly trimmed)
                                int mateLen = 0;
                                boolean trimmedPstI = false;
                                if (toks.length == 8) { //PE input
                                    builderR2.append(toks[4]); //mate id                                            
                                    if (TRIM_ADAPTERS) {
                                        int trimFrom = toks[5].indexOf("CTGCA" + SequenceOps.getReverseComplementString(barcode));
                                        if (trimFrom >= 0) {
                                            toks[5] = toks[5].substring(0, trimFrom + 5); //TODO don't trim PstI site!
                                            toks[7] = toks[7].substring(0, trimFrom + 5);
                                            trimmedPstI = true;
                                        }
                                    }
                                    builderR2.append("\t").append(toks[5]); //mate seq
                                    builderR2.append("\t").append("+"); //redundant id or '+
                                    builderR2.append("\t").append(toks[7]); //qual line
                                    mateLen = toks[5].length();
                                }
                                //If all len cutoffs met (assuming PE)
                                if (mateLen >= MIN_LENGTH_PAIR_EACH && toks[1].length() >= MIN_LENGTH_PAIR_EACH && mateLen + toks[1].length() >= MIN_LENGTH_PAIR_SUM) {
                                    bufferList.add(builderR1.append("\t").append(builderR2).toString());
                                    //else if PE input
                                } else if (toks.length == 8) {
                                    //count pairs under combined length 
                                    if (mateLen + toks[1].length() < MIN_LENGTH_PAIR_SUM) {
                                        pairUnderLenSum++;
                                    }
                                    if (toks[1].length() < MIN_LENGTH_PAIR_EACH) {
                                        pairedReadUnderLen++;
                                    } else {
                                        bufferList.add(builderR1.toString());
                                    }
                                    if (mateLen < MIN_LENGTH_PAIR_EACH) {
                                        pairedReadUnderLen++;
                                    } else {
                                        bufferList.add(builderR2.toString());
                                    }
                                    //else if SE input
                                } else if (toks.length == 4) {
                                    if (toks[1].length() >= MIN_LENGTH_READ) {
                                        bufferList.add(builderR1.toString());
                                    } else {
                                        singleUnderLength++;
                                    }
                                } else {
                                    Reporter.report("[ERROR]","FASTQ erecord expected to have 4 or 8 fileds, observed fields="+toks.length, TOOL_NAME);
                                }
                                //Cummulative trimming stats
                                if (trimmedMspI && trimmedPstI) {
                                    trimmedBothCutSitesInPair++;
                                } else if (trimmedMspI) {
                                    trimmedMspIcount++;
                                } else if (trimmedPstI) {
                                    trimmedPstIcount++;
                                }
                                break;
                            }
                        }
                        if (matchingBarcodes == 0) {
                            noBarcodeMatch++;
//                            System.out.println(line);
//                        } else if (matchingBarcodes > 1) {
//                            Reporter.report("[ERROR]", "Record matching more than one barcode " + spliPattern.split(line)[1], TOOL_NAME);
//                            System.exit(1);
                        }
                    }
                }
            }
            //output remianing stored records
            Set<String> keySet = sampleToQueueMap.keySet();

            for (String sample : keySet) {
                BlockingQueue<SampleBuffer> outQ = sampleToQueueMap.get(sample);
                SampleBuffer bufferList = sampleToBufferMap.get(sample);
                if (bufferList != null && !bufferList.isEmpty()) {
                    if (!ONLY_COUNT) {
                        outQ.put(bufferList);
                    }
                    keyMap.addToSampleCount(sample, bufferList.size());
                }
                outQ.put(new SampleBuffer());//inform other threads
            }

//            if (PRODUCER_THREADS.decrementAndGet() == 0) {
//                for (String sample : keySet) {
//                    BlockingQueue<ArrayList<String>> outQ = sampleToQueueMap.get(sample);
//                    outQ.put(new ArrayList<String>());//inform other threads
//                }
//            }
            inputQueue.put(new ArrayList<String>()); //inform other threads
            if (lines > 0) {
                String name = Thread.currentThread().getName();
                String message = "[" + name + "] " + NumberFormat.getNumberInstance().format(lines)
                    + " records processed, no matching barcode in " + NumberFormat.getNumberInstance().format(noBarcodeMatch);
                if (PSTI_STARTING_ONLY) {
                    message += ", non-PstI start detected in " + NumberFormat.getNumberInstance().format(notPstIStart);
                }
                finalMessages.add(new Message(Message.Level.INFO, message, TOOL_NAME));
                if (TRIM_ADAPTERS) {
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Trimmed total " + NumberFormat.getNumberInstance().format(trimmedMspIcount + trimmedBothCutSitesInPair + trimmedPstIcount) + " fragments", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Read-through detection with barcode and adapter trimming: ", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified MspI+3' Adapter in R1 and PstI+barcode in R2 in " + NumberFormat.getNumberInstance().format(trimmedBothCutSitesInPair) + " pairs", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified PstI+barcode in " + NumberFormat.getNumberInstance().format(trimmedPstIcount) + " R2 reads", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  Identified MspI+3' Adapter in " + NumberFormat.getNumberInstance().format(trimmedMspIcount) + " R1 reads", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "] Length filtering breakdown: ", TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(pairUnderLenSum) + " pairs under combined length " + MIN_LENGTH_PAIR_SUM, TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(pairedReadUnderLen) + " paired reads under length " + MIN_LENGTH_PAIR_EACH, TOOL_NAME));
                    finalMessages.add(new Message(Message.Level.INFO, "[" + name + "]  " + NumberFormat.getNumberInstance().format(singleUnderLength) + " single reads under length " + MIN_LENGTH_READ, TOOL_NAME));
                }
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
