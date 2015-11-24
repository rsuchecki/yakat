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
import java.util.concurrent.atomic.AtomicInteger;
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
    private final int BUFFER_SIZE = 8192; // //THAT MANY RECORDS 
    private final KeyMap keyMap;
    private final String TOOL_NAME;
    private final boolean TRIM_BARCODE;
//    private AtomicInteger PRODUCER_THREADS;
    private Integer MIN_LENGTH_READ;
    private Integer MIN_LENGTH_PAIR;

    public SplitterConsumerProducer(BlockingQueue<ArrayList<String>> inputQueue,
            KeyMap keyMap, String toolName, boolean trimBarcode,
            Integer MIN_LENGTH_READ, Integer MIN_LENGTH_PAIR) {
        this.inputQueue = inputQueue;
        this.keyMap = keyMap;
        this.TOOL_NAME = toolName;
        this.TRIM_BARCODE = trimBarcode;
//        this.PRODUCER_THREADS = producerThreads;
        this.MIN_LENGTH_READ = MIN_LENGTH_READ;
        this.MIN_LENGTH_PAIR = MIN_LENGTH_PAIR;
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

//            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    lines++;
                    tokenizer = new StringTokenizer(line);
                    String id = tokenizer.nextToken();
                    String flowcell = id.substring(1).replaceAll(":.*", "");
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
                                if (bufferList == null) { //nothing yet for this sample
                                    bufferList = new ArrayList<>(BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);
                                }
                                if (bufferList.size() == BUFFER_SIZE) { //buffer full, place on queue to write and start a new buffer
                                    BlockingQueue<ArrayList<String>> outputQueue = sampleToQueueMap.get(sample);
                                    outputQueue.put(bufferList);
                                    bufferList = new ArrayList<>(BUFFER_SIZE);
                                    sampleToBufferMap.put(sample, bufferList);
                                }
                                StringBuilder sb1 = new StringBuilder();
                                StringBuilder sb2 = new StringBuilder();
                                sb1.append(id).append("\t"); //id
                                String trimmed = sequence;
                                if (TRIM_BARCODE) {
                                    trimmed = sequence.substring(barcode.length());
                                }
                                sb1.append(trimmed); //sequence (possibly trimmed)
                                sb1.append("\t").append(tokenizer.nextToken()); //redundant id or '+
                                sb1.append("\t").append(tokenizer.nextToken()); //qual line
                                int mateLen = 0;
                                if (tokenizer.hasMoreTokens()) {
                                    sb2.append(tokenizer.nextToken()); //mate id                                            
                                    String mateSequence = tokenizer.nextToken();
                                    mateLen = mateSequence.length();
                                    sb2.append("\t").append(mateSequence);
                                    sb2.append("\t").append(tokenizer.nextToken()); //redundant id or '+
                                    sb2.append("\t").append(tokenizer.nextToken()); //qual line
                                }
                                if (mateLen >= MIN_LENGTH_READ && trimmed.length() >= MIN_LENGTH_READ && mateLen + trimmed.length() >= MIN_LENGTH_PAIR) {
                                    bufferList.add(sb1.append("\t").append(sb2).toString());
                                } else {
                                    if (trimmed.length() >= MIN_LENGTH_READ) {
                                        bufferList.add(sb1.toString());
                                    }
                                    if (mateLen >= MIN_LENGTH_READ) {
                                        bufferList.add(sb2.toString());
                                    }
                                }
                            }
                        }
                        if (matchingBarcodes == 0) {
                            noBarcodeMatch++;
//                            System.out.println(line);
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

//            if (PRODUCER_THREADS.decrementAndGet() == 0) {
//                for (String sample : keySet) {
//                    BlockingQueue<ArrayList<String>> outQ = sampleToQueueMap.get(sample);
//                    outQ.put(new ArrayList<String>());//inform other threads
//                }
//            }
            inputQueue.put(new ArrayList<String>()); //inform other threads
            if (lines > 0) {
                Reporter.report("[INFO]", "[" + Thread.currentThread().getName() + "] " + NumberFormat.getNumberInstance().format(lines) + " records processed, no matching barcode in " + NumberFormat.getNumberInstance().format(noBarcodeMatch), TOOL_NAME);
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
