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
import java.util.concurrent.BlockingQueue;
/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class WriterConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final int BUFFER_SIZE = 8192; // //THAT MANY KMERS 
    private String RECORD_NAME = "k-mers";
    private final String TOOL_NAME;

    public WriterConsumer(BlockingQueue<ArrayList<String>> outputQueue, String toolName) {
        this.outputQueue = outputQueue;
        TOOL_NAME = toolName;
    }

    public WriterConsumer(BlockingQueue<ArrayList<String>> outputQueue, String recordName, String toolName) {
        this.outputQueue = outputQueue;
        this.RECORD_NAME = recordName;
        TOOL_NAME = toolName;
    }

    @Override
    public void run() {
        BufferedWriter out = null;
        try {
//            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(java.io.FileDescriptor.out), "ASCII"), 512);
            out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(java.io.FileDescriptor.out)), BUFFER_SIZE);
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
            Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(outputCount) + " " + RECORD_NAME + " written-out",TOOL_NAME);
        } catch (InterruptedException | UnsupportedEncodingException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        } catch (IOException e) {
            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
        } finally {
            if(out != null) {
                try {
                    out.close();
                } catch (IOException ex) {
                    Reporter.report("[ERROR]", ex.getMessage(), TOOL_NAME);
                }
            }
        }
    }

}
