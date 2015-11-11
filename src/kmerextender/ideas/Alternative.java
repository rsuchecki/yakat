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
package kmerextender.ideas;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.InputReaderProducer;
import kmerextender.KmerExtender;
import kmerextender.PairMerMapPopulatorConsumer;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Alternative {

    private int KMER_LENGTH;

    public Alternative(int threads, ArrayList<String> inputFileNamesList, Integer k) {
        readKmersAndPopulateKmersMap(threads, inputFileNamesList, k);
    }

    
    
    private void readKmersAndPopulateKmersMap(int threads, ArrayList<String> inputFileNamesList, Integer k) {
//        try {
//            //READ INPUT AND POPULATE PairMers MAP
////            int threads = Math.max(Runtime.getRuntime().availableProcessors(), 6);        
//            BlockingQueue inputQueue = new ArrayBlockingQueue(65536);
////            boolean stranded = false;
//            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
//            final ExecutorService readAndPopulateDbExecutor = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
//
//            //SPAWN INPUT READING THREAD
//            InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, inputFileNamesList, k, null);
//            Future<?> future = readAndPopulateDbExecutor.submit(inputReaderProducer);
//            futures.add(future);
//            boolean splitInputSequenceIntoKmers = true;
////            if(KMER_LENGTH == null) {
////               splitInputSequenceIntoKmers = false;
////            }
//
//            //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
//            long timeStart = System.currentTimeMillis();
//            int count = 0;
//            while (inputReaderProducer.getGuessedInputFormat() == null) {
//                try {
//                    //IF nothing happens after 5 seconds
//                    if (System.currentTimeMillis() - timeStart > 5000 && (count++ % 50 == 0)) {
//                        Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...");
//                    }
//                    Thread.sleep(100); //wait for 1/10 of a second
//                } catch (InterruptedException ex) {
//                }
//            }
//            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.GuessedInputFormat.KMERS)) {
//                splitInputSequenceIntoKmers = false;
//                KMER_LENGTH = inputReaderProducer.getKmerLength();
//            }
//
//            //SPAWN THREADS TO POPULATE CLIPMERS MAP
//            for (int i = 0; i < threads; i++) {
////                KmersPopulatorConsumer consumer = new KmersPopulatorConsumer(inputQueue, splitInputSequenceIntoKmers, KMER_LENGTH, );
////                futures.add(readAndPopulateDbExecutor.submit(consumer));
//            }
//            readAndPopulateDbExecutor.shutdown();
//            try {
//                for (Future<?> f : futures) {
//                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//                }
//                readAndPopulateDbExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            } catch (InterruptedException e) {
//                Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!");
//            } catch (ExecutionException ex) {
//                Reporter.report("[ERROR]", "PairMerSet populator execution exception!");
//                ex.getMessage();
//
//            } catch (TimeoutException ex) {
//                Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
//                Reporter.report("[ERROR]", "PairMerSet populator timeout exception!");
//            }
//
////            if (pairMersMap.isOutOfMemory()) {
////                Reporter.report("[ERROR]", "Terminating, out of memory while populating the Map, k=" + KMER_LENGTH + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()));
////                System.exit(1);
////            }
//
//        } catch (OutOfMemoryError e) {
//            Reporter.report("[ERROR]", "Out of memory error!");
////            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
//            System.exit(1);
//        }

    }
}
