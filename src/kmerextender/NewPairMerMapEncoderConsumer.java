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
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class NewPairMerMapEncoderConsumer implements Runnable {

    private final BlockingQueue<List<String>> queue;
    private final NewPairMerArrays pairMerArrays;
    private final Integer MIN_KMER_FREQUENCY;
    private final boolean SPLIT_INPUT_INTO_KMERS;
//    private final Integer KMER_LENGTH;
    private final ArrayList<Integer> kList;
    private final boolean SKIP_TERMINAL_BASES;
    private final int prefixLen;
    private final int postfixLen;
    private final AtomicLongArray stats;
//    private long pairMersGenerated;

    public NewPairMerMapEncoderConsumer(BlockingQueue<List<String>> queue, NewPairMerArrays pairMerArrays,
            boolean splitInputSequenceintoKmers, ArrayList<Integer> kList, Integer minFreq, boolean skipTerminalBases, AtomicLongArray stats,
            int prefixLen, int postfixLen) {
        this.queue = queue;
        this.pairMerArrays = pairMerArrays;
        this.SPLIT_INPUT_INTO_KMERS = splitInputSequenceintoKmers;
//        this.KMER_LENGTH = k;
        this.kList = kList;
        this.MIN_KMER_FREQUENCY = minFreq;
        this.SKIP_TERMINAL_BASES = skipTerminalBases;
        this.stats = stats;
        this.prefixLen = prefixLen;
        this.postfixLen = postfixLen;
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
//                        pairMerArrays.addKmer(read, stats);
                        addKmerToMap(read, 1);
//                        if (!addKmerToMap(read, 1)) {
//                            break;
//                        }
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
//                        pairMerArrays.addKmer(read, stats);
                        addKmerToMap(read, 1);
//                        if (!addKmerToMap(read, frequency)) {
//                            break;
//                        }
                    }
                }
            }
            queue.put(new ArrayList<String>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
//        pairMersMaps.addToTotalPairMersGenerated(pairMersGenerated);
    }

    private boolean addKmerToMap(CharSequence inputKmerSequence, int freq) throws InterruptedException {
        try {
//            System.err.println("");
//            System.err.println(" "+kmerString+" <-Generating from");
//            synchronized (KmerExtender.class) {
//                pairMerArrays.getPairMersMap(kmerString.length()).addToPairMersMap(kmerString, 0, kmerString.length() - 1, true, !SPLIT_INPUT_INTO_KMERS, freq);
//                pairMerArrays.getPairMersMap(kmerString.length()).addToPairMersMap(kmerString, 0, kmerString.length() - 1, false, !SPLIT_INPUT_INTO_KMERS, freq);
//                  pairMerArrays.addKmer(kmerString, stats);

            //Process each of the 2  paiMers separately
            //Separate left clip+core            
            char leftClip = inputKmerSequence.charAt(0);
            CharSequence core1 = inputKmerSequence.subSequence(1, inputKmerSequence.length());
            boolean is1rc = false;
            if (!SequenceOps.isCanonical(core1)) {
                core1 = SequenceOps.getReverseComplement(core1);
                leftClip = SequenceOps.complement(leftClip);
                is1rc = true;
            }
                //        int prefix = CoreCoder.encodeCore(core.subSequence(0, prefixLen));
                pairMerArrays.addMer(core1.toString(), leftClip, !is1rc, stats);

            //Separate core+right clip
            CharSequence core2 = inputKmerSequence.subSequence(0, inputKmerSequence.length() - 1);
            char rightClip = inputKmerSequence.charAt(inputKmerSequence.length() - 1);
            boolean is2rc = false;
            if (!SequenceOps.isCanonical(core2)) {
                core2 = SequenceOps.getReverseComplement(core2);
                rightClip = SequenceOps.complement(rightClip);
                is2rc = true;
            }

                pairMerArrays.addMer(core2.toString(), rightClip, is2rc, stats);

//            }
//            pairMersGenerated+=2;
        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error !", this.getClass().getSimpleName());
//            pairMerArrays.getPairMersMap(kmerString.length()).setOutOfMemory();
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

    private void kmerizeAndAddToMaps(CharSequence sequence) throws InterruptedException {
        for (Integer k : kList) {
            try {
                int maxKmer = sequence.length() - k + 1;
                int startAt = 0;
                if (SKIP_TERMINAL_BASES) { //to exclude terminal k-mers
                    --maxKmer;
                    startAt++;
                }
                for (int i = startAt; i < maxKmer; i++) {
//                String s = (String) sequence.subSequence(i, i + KMER_LENGTH);
//                    pairMerArrays.addKmer(sequence.subSequence(i, i + k), stats);
                    addKmerToMap(sequence.subSequence(i, i + k), 1);

//                    pairMerArrays.getPairMersMap(k).addToPairMersMap(sequence, i, i + k - 1, true, !SPLIT_INPUT_INTO_KMERS, 1);
//                    pairMerArrays.getPairMersMap(k).addToPairMersMap(sequence, i, i + k - 1, false, !SPLIT_INPUT_INTO_KMERS, 1);
                }
            } catch (OutOfMemoryError e) {
                Reporter.report("[ERROR]", "Out of memory error !", this.getClass().getSimpleName());
//                pairMerArrays.getPairMersMap(k).setOutOfMemory();
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
