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
package processpileup;

import gbssplit.*;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PileupConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> inputQueue;
    private final int OUT_BUFFER_SIZE;
    private final String TOOL_NAME;

    public PileupConsumer(BlockingQueue<ArrayList<String>> inputQueue,
             String toolName, int OUT_BUFFER_SIZE) {
        this.inputQueue = inputQueue;
        this.TOOL_NAME = toolName;
//        this.PRODUCER_THREADS = producerThreads;
        this.OUT_BUFFER_SIZE = OUT_BUFFER_SIZE;
    }

    @Override
    public void run() {
        try {
            ArrayList<String> list;
            long noBarcodeMatch = 0L;
            long lines = 0L;
            while (!(list = inputQueue.take()).isEmpty()) {
                for (String line : list) {
                    lines++;
                    StringTokenizer tokenizer = new StringTokenizer(line);
                    String id = tokenizer.nextToken();
                    int position = Integer.parseInt(tokenizer.nextToken());
                    char refBase = tokenizer.nextToken().charAt(0);
                    while(tokenizer.hasMoreTokens()) {
                        int coverage = Integer.parseInt(tokenizer.nextToken());
//                        tokenizer.nextToken();
//                        tokenizer.nextToken();
                        String bases = tokenizer.nextToken();
                        String quals = tokenizer.nextToken();
                    }
                }
            }
            //output remianing stored records
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
