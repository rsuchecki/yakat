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
package kmermatch;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import shared.InputReaderProducer.InFormat;
import shared.Message;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerMatcherConsumerProducer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final ConcurrentSkipListSet<Kmer> map;
    private final boolean invertMatch;
    private final String TOOL_NAME;
    private final ArrayList<Message> finalMessages;
    private final int BUFFER_SIZE;
    private final InFormat inFormat;
    private final int k;

    public KmerMatcherConsumerProducer(BlockingQueue<ArrayList<String>> inputQueue, BlockingQueue<ArrayList<String>> outputQueue,
            ConcurrentSkipListSet<Kmer> kmers, int BUFFER_SIZE, String TOOL_NAME, ArrayList<Message> finalMessages, boolean invertMatch, InFormat inFormat, int k) {
        this.inputQueue = inputQueue;
        this.outputQueue = outputQueue;
        this.map = kmers;
        this.TOOL_NAME = TOOL_NAME;
        this.finalMessages = finalMessages;
        this.BUFFER_SIZE = BUFFER_SIZE;
        this.invertMatch = invertMatch;
        this.inFormat = inFormat;
        this.k = k;
    }

    @Override
    public void run() {
        try {
            Pattern spliPattern = Pattern.compile("\t");
            ArrayList<String> list;
            ArrayList<String> buffer = new ArrayList<>(BUFFER_SIZE);

            long countIn = 0;
            long countOut = 0;
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    String toks[] = spliPattern.split(line);
                    //Currently only considering one record per line wrapped FAST(A|Q) PE or SE
//                    if (toks.length == 2 || (toks.length == 4 && inFormat == InFormat.FASTQ_SE_ONE_LINE)) {  //FAST(A|Q) SE 
//                        //kmerize toks[1]                        
//                    } else if (toks.length == 8) { //FASTQ PE
//                        //kmerize toks[1] and toks[4]                        
//                    } else { //FASTA PE
//                        //kmerize toks[1] and toks[3]
//                    }

                    
                    
                    int minMatches = 1;
                    int matches = kmerMatches(toks[1]);
                    if (matches >= minMatches && !invertMatch || matches < minMatches && invertMatch) {
                        if (buffer.size() >= BUFFER_SIZE) {
                            putOneQueue(outputQueue, buffer);
                            buffer = new ArrayList<>();
                        }
                        countIn++;
                        buffer.add(line);
                    } else {
                        countOut++;
                    }
                }
            }
            if (!buffer.isEmpty()) {
                putOneQueue(outputQueue, buffer);
            }
            putOneQueue(inputQueue, new ArrayList<>(0)); //inform other threads                
            putOneQueue(outputQueue, new ArrayList<>(0)); //inform other threads                
            System.err.println(countIn + "=In  Out=" + countOut);
        } catch (InterruptedException ex) {
            Logger.getLogger(KmerMatcherConsumerProducer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

//    private boolean matchKmers(int minMatches, CharSequence sequence) {
//        return kmerMatches(sequence) >= minMatches;
//    }
//    
//    private boolean matchKmers(int minMatches, CharSequence sequence1, CharSequence sequence2) {       
//        return kmerMatches(sequence1)+kmerMatches(sequence2) >= minMatches;
//    }
    
    private int kmerMatches(CharSequence sequence) {
        int startAt = 0;
        int maxKmer = sequence.length() - k + 1;
        int matches = 0;
//        System.err.println(sequence);
//        String format = "%s\n";
        for (int i = startAt; i < maxKmer; i++) {
//            if(i>0) {
//                format = "%"+i+"s%s\n";
//                System.err.printf(format, "",sequence.subSequence(i, i + k));            
//            } else
//                System.err.printf(format, sequence.subSequence(i, i + k));            
            if(map.contains(new KmerASCII(SequenceOps.getCanonical(sequence.subSequence(i, i+k).toString())))) {
                matches++;
            } 
        }
//        System.exit(1);
        return matches;
    }

    private void putOneQueue(BlockingQueue q, List buffer) throws InterruptedException {
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode());
        q.put(buffer);
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode() + " done");
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode() + " first record " + buffer.get(0));
    }
}
