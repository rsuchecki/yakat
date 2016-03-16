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
package kmermatch.TODO;

import kmerextender.*;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerMapPopulatorConsumer implements Runnable {

    private final BlockingQueue<List<String>> queue;
    private final KmerMap kmerMap;
    private final Integer MIN_KMER_FREQUENCY;
    private final boolean SPLIT_INPUT_INTO_KMERS;
    private final String TOOL_NAME;
    private final Integer KMER_LENGTH;

    public KmerMapPopulatorConsumer(BlockingQueue<List<String>> queue, KmerMap kmerMap,
        boolean splitInputSequenceintoKmers, Integer minFreq, String TOOL_NAME, Integer k) {
        this.queue = queue;
        this.kmerMap = kmerMap;
        this.SPLIT_INPUT_INTO_KMERS = splitInputSequenceintoKmers;
//        this.KMER_LENGTH = k;
        this.MIN_KMER_FREQUENCY = minFreq;
        this.TOOL_NAME = TOOL_NAME;
        this.KMER_LENGTH = k;
    }

    @Override
    public void run() {
        try {
            StringTokenizer tokenizer;
            List<String> list;
            if (SPLIT_INPUT_INTO_KMERS) {
                while (!(list = queue.take()).isEmpty()) {
                    for (String read : list) {
                        kmerizeAndAddToMap(read);
                    }
                }
            } else { //SIMPLY LOADING A BUNCH OF KMERS 
                while (!(list = queue.take()).isEmpty()) {
                    for (String read : list) {
                        tokenizer = new StringTokenizer(read);
                        read = tokenizer.nextToken().toUpperCase();
                        if (MIN_KMER_FREQUENCY != null && (tokenizer.hasMoreTokens())) {
                            int frequency = Integer.parseInt(tokenizer.nextToken());
                            if (frequency < MIN_KMER_FREQUENCY) {
                                continue;
                            }
                        }
                        if (!addKmerToMap(read)) {
                            break;
                        }
                    }
                }
            }
            queue.put(new ArrayList<String>()); //inform other threads

        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private boolean addKmerToMap(CharSequence kmerString) {
        try {
//            System.err.println("");
//            System.err.println(" "+kmerString+" <-Generating from");
//            pairMersMap.addToPairMersMap(kmerString, true, overlapLength, !SPLIT_INPUT_INTO_KMERS);
//            pairMersMap.addToPairMersMap(kmerString, false, overlapLength, !SPLIT_INPUT_INTO_KMERS);
            kmerMap.addToKmersMap(kmerString, 0, kmerString.length() - 1, true);

        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error !", TOOL_NAME);
//            kmerMap.getPairMersMap(kmerString.length()).setOutOfMemory();
//            try {
//                queue.put(new ArrayList<String>());
//            } catch (InterruptedException ex) {
//                System.err.println(ex.getMessage());
//            }
////            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            return false;
        }
        return true;
    }

    private void kmerizeAndAddToMap(CharSequence sequence) {
            try {
                int maxKmer = sequence.length() - KMER_LENGTH + 1;
                for (int i = 0; i < maxKmer; i++) {
//                String s = (String) sequence.subSequence(i, i + KMER_LENGTH);
                    kmerMap.addToKmersMap(sequence, i, i + KMER_LENGTH, true);
                }
            } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error !", TOOL_NAME);
//                kmerMap.getPairMersMap(k).setOutOfMemory();
//                try {
//                    queue.put(new ArrayList<String>());
//                } catch (InterruptedException ex) {
//                    System.err.println(ex.getMessage());
//                }
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            }
    }

}
