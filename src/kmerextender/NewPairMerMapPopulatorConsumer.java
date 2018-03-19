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
import shared.Reporter;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class NewPairMerMapPopulatorConsumer implements Runnable {

    private final NewPairMerArrays pairMerArrays;
    private final List<BlockingQueue<List<NewPairMerWithPrefixes>>> queues;
    
    public NewPairMerMapPopulatorConsumer(NewPairMerArrays pairMerArrays, List<BlockingQueue<List<NewPairMerWithPrefixes>>> queues) {
        this.pairMerArrays = pairMerArrays;
        this.queues = queues;
    }
    
    

    @Override
    public void run() {

//        try {
//            
//            for (int i = 0; i < queues.size(); i++) {
//                BlockingQueue<List<NewPairMerWithPrefixes>> queue = queues.get(i);
//                List<NewPairMerWithPrefixes> list;                
//                while (!(list = queue.take()).isEmpty()) {
//                    for (NewPairMerWithPrefixes newPairMerWithPrefixes : list) {
//                        
//                    }
//                }
//                
//            }
////            List<String> list;
//            
//            for (int i = 0; i < queues.size(); i++) {
//                queue.put(new ArrayList<String>()); //inform other threads
//            }
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }
    }
}