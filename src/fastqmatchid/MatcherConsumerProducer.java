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
package fastqmatchid;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import shared.Message;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class MatcherConsumerProducer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final ConcurrentSkipListSet<Identifier> ids;
    private final boolean invertMatch;
    private final String TOOL_NAME;
    private final ArrayList<Message> finalMessages;
    private final int BUFFER_SIZE;

    public MatcherConsumerProducer(BlockingQueue<ArrayList<String>> inputQueue, BlockingQueue<ArrayList<String>> outputQueue,
            ConcurrentSkipListSet<Identifier> ids, int BUFFER_SIZE, String TOOL_NAME, ArrayList<Message> finalMessages, boolean invertMatch) {
        this.inputQueue = inputQueue;
        this.outputQueue = outputQueue;
        this.ids = ids;
        this.TOOL_NAME = TOOL_NAME;
        this.finalMessages = finalMessages;
        this.BUFFER_SIZE = BUFFER_SIZE;
        this.invertMatch = invertMatch;
    }

    @Override
    public void run() {
        try {
            Pattern spliPattern = Pattern.compile("\t| |/");
            ArrayList<String> list;
            ArrayList<String> buffer = new ArrayList<>(BUFFER_SIZE);

//            long countIn = 0;
//            long countOut = 0;
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    String toks[] = spliPattern.split(line);
//                    countIn++;
//                    String replaceFirst = toks[0].replaceFirst("\\/[1-2]", "");
                    
                    boolean contains = ids.contains(new Identifier(toks[0].substring(1)));//.replaceFirst("\\/[1-2]", "")));
                    if (contains && !invertMatch || !contains && invertMatch) {
                        if (buffer.size() >= BUFFER_SIZE) {
                            putOneQueue(outputQueue, buffer);
                            buffer = new ArrayList<>();
                        }
                        buffer.add(line);
                    }
                }
            }
            if (!buffer.isEmpty()) {
                putOneQueue(outputQueue, buffer);
            }
            putOneQueue(inputQueue, new ArrayList<>(0)); //inform other threads                
            putOneQueue(outputQueue, new ArrayList<>(0)); //inform other threads                
//            System.err.println(countIn + "=In  Out=" + countOut);
        } catch (InterruptedException ex) {
            Logger.getLogger(MatcherConsumerProducer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private void putOneQueue(BlockingQueue q, List buffer) throws InterruptedException {
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode());
        q.put(buffer);
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode() + " done");
//        System.err.println(Thread.currentThread().getName() + "[MC] puts " + buffer.size() + " on " + q.hashCode() + " first record " + buffer.get(0));
    }
}
