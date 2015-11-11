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
package shared;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.Reporter;
/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class WriterConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final int BUFFER_SIZE = 8192; // //THAT MANY KMERS 
    private String recordName = "k-mers";

    public WriterConsumer(BlockingQueue<ArrayList<String>> outputQueue) {
        this.outputQueue = outputQueue;
    }

    public WriterConsumer(BlockingQueue<ArrayList<String>> outputQueue, String recordName) {
        this.outputQueue = outputQueue;
        this.recordName = recordName;
    }

    @Override
    public void run() {
        try {
//            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(java.io.FileDescriptor.out), "ASCII"), 512);
            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(java.io.FileDescriptor.out)), BUFFER_SIZE);
            ArrayList<String> list;
            long outputCount = 0L;
            while (!(list = outputQueue.take()).isEmpty()) {
                outputCount += list.size();
                for (String line : list) {
                    out.write(line);
                    out.newLine();
                }
            }
            out.flush();
            //output remianing stored records
            outputQueue.put(new ArrayList<String>());            
            Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(outputCount) + " " + recordName + " written-out", getClass().getSimpleName());
        } catch (InterruptedException | UnsupportedEncodingException e) {
            Reporter.report("[ERROR]", e.getMessage(), getClass().getSimpleName());
        } catch (IOException e) {
            Reporter.report("[ERROR]", e.getMessage(), getClass().getSimpleName());
        }
    }

}
