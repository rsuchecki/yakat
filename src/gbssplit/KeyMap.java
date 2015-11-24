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
package gbssplit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.GZIPInputStream;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KeyMap {

    private ConcurrentHashMap<String, ConcurrentHashMap<String, String>> keyMap;
    private ConcurrentHashMap<String, BlockingQueue<ArrayList<String>>> sampleToQueueMap;
    private final int READER_BUFFER_SIZE = 8192;
    private final int QUEUE_BUFFER_SIZE = 255;
    private final String TOOL_NAME;
    private final String BLANK_SAMPLE_NAME;

    public KeyMap(String keyFileName, String toolName, String BlankSampleName) {
        TOOL_NAME = toolName;
        populateMap(keyFileName);
        BLANK_SAMPLE_NAME = BlankSampleName;
    }

    public String getSample(String flowCell, String barcode) {
        return keyMap.get(flowCell).get(barcode);
    }
    
    public ConcurrentHashMap getBarcodes(String flowCell) {
        return keyMap.get(flowCell);
    }
        
    public int getSamplesTotal() {
        int tot = 0;
        for (Map.Entry<String, ConcurrentHashMap<String, String>> entrySet : keyMap.entrySet()) {            
            tot += entrySet.getValue().size();
        }
        return tot;
    }
    
    private void populateMap(String inputFile) {
        keyMap = new ConcurrentHashMap<>(1000);
        sampleToQueueMap = new ConcurrentHashMap<>(1000);
        BufferedReader content = null;
        try {
            if (inputFile.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile), READER_BUFFER_SIZE);
                content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                content = new BufferedReader(new FileReader(new File(inputFile)), READER_BUFFER_SIZE);
            }
            String line;
            while ((line = content.readLine()) != null && !line.isEmpty()) {
                String[] toks = line.split("\t");
                String flowcell = toks[0];
                String barcode = toks[2];
                String sample = toks[3];
                if(sample.equalsIgnoreCase(BLANK_SAMPLE_NAME)) {
                    for (int i = 4; i < toks.length; i++) {
                        sample += "_"+toks[i];
                    }
                }
                if (keyMap.containsKey(flowcell)) {
                    ConcurrentHashMap<String, String> barcodeToSample = keyMap.get(flowcell);
                    if (barcodeToSample.containsKey(barcode)) {
                        Reporter.report("[ERROR]", "Duplicate barcode on a flowcell!!!", TOOL_NAME);
                        System.exit(1);
                    } else {
                        barcodeToSample.put(barcode, sample);
                        sampleToQueueMap.put(sample, new ArrayBlockingQueue<ArrayList<String>>(QUEUE_BUFFER_SIZE));
                    }
                } else {
                    ConcurrentHashMap<String, String> barcodeToSample = new ConcurrentHashMap<>();
                    barcodeToSample.put(barcode, sample);
                    keyMap.put(flowcell, barcodeToSample);
                    sampleToQueueMap.put(sample, new ArrayBlockingQueue<ArrayList<String>>(QUEUE_BUFFER_SIZE));
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found exception: " + ex.getMessage(), TOOL_NAME);
            System.exit(1);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        } finally {
            try {
                if (content != null) {
                    content.close();
                }
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
    }

    public ConcurrentHashMap<String, BlockingQueue<ArrayList<String>>> getSampleToQueueMap() {
        return sampleToQueueMap;
    }

    
}
