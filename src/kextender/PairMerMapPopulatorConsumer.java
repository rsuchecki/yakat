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
package kextender;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class PairMerMapPopulatorConsumer implements Runnable {

    private final BlockingQueue<List<String>> queue;
    private final PairMerMaps pairMersMaps;
    private final int MIN_KMER_FREQUENCY;
    private final boolean SPLIT_INPUT_INTO_KMERS;
//    private final Integer KMER_LENGTH;
    private final ArrayList<Integer> kList;
    private final boolean SKIP_TERMINAL_BASES;
//    private long pairMersGenerated;

    public PairMerMapPopulatorConsumer(BlockingQueue<List<String>> queue, PairMerMaps pairMerMaps,
            boolean splitInputSequenceintoKmers, ArrayList<Integer> kList, int minFreq, boolean skipTerminalBases) {
        this.queue = queue;
        this.pairMersMaps = pairMerMaps;
        this.SPLIT_INPUT_INTO_KMERS = splitInputSequenceintoKmers;
//        this.KMER_LENGTH = k;
        this.kList = kList;
        this.MIN_KMER_FREQUENCY = minFreq;
        this.SKIP_TERMINAL_BASES = skipTerminalBases;
    }

    @Override
    public void run() {
        try {
            StringTokenizer tokenizer;
            List<String> list;
            if (SPLIT_INPUT_INTO_KMERS) {
                while (!(list = queue.take()).isEmpty()) {
                    for (String read : list) {
                        kmerizeAndAddToMaps(read);
                    }
                }
            } else if (MIN_KMER_FREQUENCY <= 1) {
                while (!(list = queue.take()).isEmpty()) {
                    for (String read : list) {
                        tokenizer = new StringTokenizer(read);
                        read = tokenizer.nextToken().toUpperCase();
                        if (!addKmerToMap(read, 1)) {
                            break;
                        }
                    }
                }
            } else {
                while (!(list = queue.take()).isEmpty()) {
                    for (String read : list) {
                        tokenizer = new StringTokenizer(read);
                        read = tokenizer.nextToken().toUpperCase();
                        int frequency = 1;
                        if (tokenizer.hasMoreTokens()) {
                            frequency = Integer.parseInt(tokenizer.nextToken());
                            if (frequency < MIN_KMER_FREQUENCY) {
                                continue;
                            }
                        }
                        if (!addKmerToMap(read, frequency)) {
                            break;
                        }
                    }
                }
            }
            queue.put(new ArrayList<String>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
//        pairMersMaps.addToTotalPairMersGenerated(pairMersGenerated);
    }

    private boolean addKmerToMap(CharSequence kmerString, int freq) {
        try {
//            System.err.println("");
//            System.err.println(" "+kmerString+" <-Generating from");
//            synchronized (KmerExtender.class) {
                pairMersMaps.getPairMersMap(kmerString.length()).addToPairMersMap(kmerString, 0, kmerString.length() - 1, true, !SPLIT_INPUT_INTO_KMERS, freq);
                pairMersMaps.getPairMersMap(kmerString.length()).addToPairMersMap(kmerString, 0, kmerString.length() - 1, false, !SPLIT_INPUT_INTO_KMERS, freq);
//            }
//            pairMersGenerated+=2;
        } catch (OutOfMemoryError e) {
//            Reporter.report("[ERROR]", "Out of memory error !");
            pairMersMaps.getPairMersMap(kmerString.length()).setOutOfMemory();
            try {
                queue.put(new ArrayList<>());
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            return false;
        }
        return true;
    }

    private void kmerizeAndAddToMaps(CharSequence sequence) {
        for (int k : kList) {
            try {
                int maxKmer = sequence.length() - k + 1;
                int startAt = 0;
                if (SKIP_TERMINAL_BASES) { //to exclude terminal k-mers
                    --maxKmer;
                    startAt++;
                }
                for (int i = startAt; i < maxKmer; i++) {
//                String s = (String) sequence.subSequence(i, i + KMER_LENGTH);
                    pairMersMaps.getPairMersMap(k).addToPairMersMap(sequence, i, i + k - 1, true, !SPLIT_INPUT_INTO_KMERS, 1);
                    pairMersMaps.getPairMersMap(k).addToPairMersMap(sequence, i, i + k - 1, false, !SPLIT_INPUT_INTO_KMERS, 1);
                }
            } catch (OutOfMemoryError e) {
//            Reporter.report("[ERROR]", "Out of memory error !");
                pairMersMaps.getPairMersMap(k).setOutOfMemory();
                try {
                    queue.put(new ArrayList<String>());
                } catch (InterruptedException ex) {
                    System.err.println(ex.getMessage());
                }
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            }
        }
    }
}
