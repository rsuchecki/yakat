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

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SplitterConsumerProducer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
//    private final Integer MIN_KMER_FREQUENCY;
//    private final Integer MAX_KMER_FREQUENCY;
//    private final String OUT_LABEL;
    private final int BUFFER_SIZE = 100; // //THAT MANY KMERS 
    private final KeyMap keyMap;
    private final String TOOL_NAME;
    private final boolean TRIM_BARCODE;

    public SplitterConsumerProducer(BlockingQueue<ArrayList<String>> inputQueue,
            KeyMap keyMap, String toolName, boolean trimBarcode) {
        this.inputQueue = inputQueue;
        this.keyMap = keyMap;
        this.TOOL_NAME = toolName;
        this.TRIM_BARCODE = trimBarcode;
    }

    @Override
    public void run() {
        try {
            ConcurrentHashMap<String, ArrayList<String>> sampleToBufferMap = new ConcurrentHashMap<>(keyMap.getSamplesTotal() * 5);
            ConcurrentHashMap<String, BlockingQueue<ArrayList<String>>> sampleToQueueMap = keyMap.getSampleToQueueMap();
            StringTokenizer tokenizer;
            ArrayList<String> list;
            long noBarcodeMatch = 0L;
            long lines = 0L;
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    lines++;
                    tokenizer = new StringTokenizer(line);
                    String flowcell = tokenizer.nextToken().substring(1).replaceAll(":.*", "");
                    ConcurrentHashMap<String, String> barcodes = keyMap.getBarcodes(flowcell);
                    if (barcodes != null) {
                        String sequence = tokenizer.nextToken();
                        int matchingBarcodes = 0;
                        for (Map.Entry<String, String> entrySet : barcodes.entrySet()) {
                            String barcode = entrySet.getKey();
                            if (sequence.startsWith(barcode)) {
                                matchingBarcodes++;
                                String sample = entrySet.getValue();
                                ArrayList<String> bufferList = sampleToBufferMap.get(sample);
                                if (bufferList == null) {
                                    bufferList = new ArrayList<>(BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);
                                }
                                if (bufferList.size() == BUFFER_SIZE) {
                                    BlockingQueue<ArrayList<String>> outputQueue = sampleToQueueMap.get(sample);
                                    outputQueue.put(bufferList);
                                    bufferList = new ArrayList<>(BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);
                                }
                                if (TRIM_BARCODE) {
                                    StringBuilder sb = new StringBuilder();
                                    tokenizer = new StringTokenizer(line);
                                    sb.append(tokenizer.nextToken()).append("\t"); //id
                                    sb.append(tokenizer.nextToken().substring(barcode.length())); //trimmed sequence
                                    while (tokenizer.hasMoreTokens()) {
                                        sb.append("\t").append(tokenizer.nextToken());
                                    }
                                    bufferList.add(line);
                                } else {
                                    bufferList.add(line);
                                }
                            }
                        }
                        if (matchingBarcodes == 0) {
                            noBarcodeMatch++;
                        } else if (matchingBarcodes > 1) {
                            Reporter.report("[ERROR]", "Sequence matching more than one barcode: " + sequence, TOOL_NAME);
                            System.exit(1);
                        }
                    }
                }
            }
            //output remianing stored records
            Set<String> keySet = sampleToQueueMap.keySet();
            for (String sample : keySet) {
                BlockingQueue<ArrayList<String>> outQ = sampleToQueueMap.get(sample);
                ArrayList<String> bufferList = sampleToBufferMap.get(sample);
                if (bufferList != null && !bufferList.isEmpty()) {
                    outQ.put(bufferList);
                }
                outQ.put(new ArrayList<String>());//inform other threads
            }
            inputQueue.put(new ArrayList<String>()); //inform other threads
            Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(lines)+" records processed, no matching barcode in "+NumberFormat.getNumberInstance().format(noBarcodeMatch), TOOL_NAME);

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
