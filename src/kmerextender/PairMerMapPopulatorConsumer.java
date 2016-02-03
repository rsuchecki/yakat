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
package kmerextender;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMerMapPopulatorConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> queue;
    private final PairMersMap pairMersMap;
    private final Integer MIN_KMER_FREQUENCY;
    private final boolean SPLIT_INPUT_INTO_KMERS;
    private final Integer KMER_LENGTH;
    private final boolean SKIP_TERMINAL_BASES;

    public PairMerMapPopulatorConsumer(BlockingQueue<ArrayList<String>> queue, PairMersMap pairMersMap,
        boolean splitInputSequenceintoKmers, Integer k, Integer minFreq, boolean skipTerminalBases) {
        this.queue = queue;
        this.pairMersMap = pairMersMap;
        this.SPLIT_INPUT_INTO_KMERS = splitInputSequenceintoKmers;
        this.KMER_LENGTH = k;
        this.MIN_KMER_FREQUENCY = minFreq;
        this.SKIP_TERMINAL_BASES = skipTerminalBases;
    }

    @Override
    public void run() {
        try {
            StringTokenizer tokenizer;
            ArrayList<String> list;
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
            pairMersMap.addToPairMersMap(kmerString, 0, kmerString.length()-1, true, !SPLIT_INPUT_INTO_KMERS);
            pairMersMap.addToPairMersMap(kmerString, 0, kmerString.length()-1, false, !SPLIT_INPUT_INTO_KMERS);

        } catch (OutOfMemoryError e) {
//            Reporter.report("[ERROR]", "Out of memory error !");
            pairMersMap.setOutOfMemory();
            try {
                queue.put(new ArrayList<String>());
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            return false;
        }
        return true;
    }

    private void kmerizeAndAddToMap(CharSequence sequence) {
        try {
            int maxKmer = sequence.length() - KMER_LENGTH + 1;
            int startAt = 0;
            if (SKIP_TERMINAL_BASES) { //to exclude terminal k-mers
                --maxKmer;
                startAt++;
            }
            for (int i = startAt; i < maxKmer; i++) {
//                String s = (String) sequence.subSequence(i, i + KMER_LENGTH);
                pairMersMap.addToPairMersMap(sequence, i, i + KMER_LENGTH - 1, true, !SPLIT_INPUT_INTO_KMERS);
                pairMersMap.addToPairMersMap(sequence, i, i + KMER_LENGTH - 1, false, !SPLIT_INPUT_INTO_KMERS);
            }
        } catch (OutOfMemoryError e) {
//            Reporter.report("[ERROR]", "Out of memory error !");
            pairMersMap.setOutOfMemory();
            try {
                queue.put(new ArrayList<String>());
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
        }
    }

//    private ArrayList<String> getKmers(String sequence) {
//        int maxKmer = sequence.length() - KMER_LENGTH + 1;
//        int startAt = 0;
//        if (SKIP_TERMINAL_BASES) { //to exclude terminal k-mers
//            --maxKmer;
//            startAt++;
//        }
//        if (maxKmer < 0) { //in case of short sequences
//            return new ArrayList<>(0);
//        }
//        ArrayList<String> kmers = new ArrayList<>(maxKmer);
//        for (int i = startAt; i < maxKmer; i++) {
//            String s = sequence.substring(i, i + KMER_LENGTH);
//            kmers.add(s);
//        }
//        return kmers;
//    }
}
