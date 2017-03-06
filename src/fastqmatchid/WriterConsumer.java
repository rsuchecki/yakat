/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.NumberFormat;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class WriterConsumer implements Runnable {

    private final BlockingQueue<List<String>> outputQueue;
    private final String outFile;
    private int matcherThreads;
    private final int WRITER_BUFFER_SIZE = 8192;
    private final String TOOL_NAME;

    public WriterConsumer(BlockingQueue<List<String>> outputQueue, String outFile, int matcherThreads, String TOOL_NAME) {
        this.outputQueue = outputQueue;
        this.outFile = outFile;
        this.matcherThreads = matcherThreads;        
        this.TOOL_NAME = TOOL_NAME;
    }
    
    

    @Override
    public void run() {
        long count = 0;
        BufferedWriter out = null;
        try {
            if (outFile.endsWith(".gz")) {
                out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile)), "UTF-8"), WRITER_BUFFER_SIZE);
            } else {
                out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"), WRITER_BUFFER_SIZE);
            }
        } catch (IOException e) {
            Logger.getLogger(WriterConsumer.class.getName()).log(Level.SEVERE, null, e);
        }

        try {
            List<String> list;
            while (!(list = outputQueue.take()).isEmpty() || --matcherThreads > 0) {
                count+=list.size();
                for (String line : list) {
                    out.write(line);
                    out.newLine();
                }
                out.flush();
            }
        } catch (InterruptedException | IOException ex) {
            Logger.getLogger(WriterConsumer.class.getName()).log(Level.SEVERE, null, ex);
        }
        Reporter.report("[INFO]", "Finished outputing matched (or unmatched) records, n=" + NumberFormat.getNumberInstance().format(count), TOOL_NAME);

    }

}
