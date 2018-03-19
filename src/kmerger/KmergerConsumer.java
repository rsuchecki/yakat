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
package kmerger;

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmergerConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final Integer MIN_KMER_FREQUENCY;
    private final Integer MAX_KMER_FREQUENCY;
    private final String OUT_LABEL;
    private final int KMER_BUFFER_SIZE = 8192; // //THAT MANY KMERS 

    public KmergerConsumer(BlockingQueue<ArrayList<String>> inputQueue, BlockingQueue<ArrayList<String>> outputQueue,
            Integer minFreq, Integer maxFreq, String outputLabel) {
        this.inputQueue = inputQueue;
        this.outputQueue = outputQueue;
        MIN_KMER_FREQUENCY = minFreq;
        MAX_KMER_FREQUENCY = maxFreq;
        OUT_LABEL = outputLabel;
    }

    @Override
    public void run() {
        try {
            ArrayList<String> bufferList = new ArrayList<>(KMER_BUFFER_SIZE);
            StringTokenizer tokenizer;
            ArrayList<String> list;
            //SIMPLY LOADING A BUNCH OF KMERS 
            String previous = null;
            Long previousFreq = null;
            String previousLabel = null;
            while (!(list = inputQueue.take()).isEmpty()) {
//                    if (!read.isEmpty()) {
                for (String line : list) {
                    tokenizer = new StringTokenizer(line);
                    String kmer = tokenizer.nextToken().toUpperCase();
                    if (previous == null) { //first k-mer
                        previous = kmer;
                        if (tokenizer.hasMoreTokens()) {
                            previousFreq = Long.parseLong(tokenizer.nextToken());
                        }
                        if (tokenizer.hasMoreTokens()) {
                            previousLabel = tokenizer.nextToken();
                        }
                    } else { //subsequent k-mers
                        if (kmer.equals(previous)) { //duplicate
                            if (tokenizer.hasMoreTokens()) {
                                previousFreq += Integer.parseInt(tokenizer.nextToken());
                            }
                            if (tokenizer.hasMoreTokens()) {
//                                String Label = tokenizer.nextToken();
                                //do we compare the labels?
                            }
                        } else { //different
                            //output previous
                            addRecordToList(previous, previousFreq, previousLabel, bufferList);                            
                            if (bufferList.size() == KMER_BUFFER_SIZE) {
                                outputQueue.put(bufferList);
                                bufferList = new ArrayList<>();
                            }
                            previous = kmer;
                            if (tokenizer.hasMoreTokens()) {
                                previousFreq = Long.parseLong(tokenizer.nextToken());
                            }
                        }

                    }
                }
            }
            //store last record 
            addRecordToList(previous, previousFreq, previousLabel, bufferList);
            //output remianing stored records
            outputQueue.put(bufferList);
            inputQueue.put(new ArrayList<String>()); //inform other threads
            outputQueue.put(new ArrayList<String>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private void addRecordToList(String previous, Long previousFreq, String previousLabel, ArrayList<String> bufferList) {
        if ((MIN_KMER_FREQUENCY == null || previousFreq >= MIN_KMER_FREQUENCY) && (MAX_KMER_FREQUENCY == null || previousFreq <= MAX_KMER_FREQUENCY)) {
            StringBuilder sb = new StringBuilder(previous);
            if (previousFreq != null) {
                sb.append("\t").append(previousFreq);
            }
            if (OUT_LABEL != null) {
                sb.append("\t").append(OUT_LABEL);
            } else if (previousLabel != null) {
                sb.append("\t").append(previousLabel);
            }
            bufferList.add(sb.toString());
        }
    }

}
