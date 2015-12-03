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
    private ConcurrentHashMap<String, BlockingQueue<SampleBuffer>> sampleToQueueMap;
    private ConcurrentHashMap<String, Long> sampleToCountMap;
    
    private final int OUT_Q_CAPACITY;
    private final String TOOL_NAME;
    private final String BLANK_SAMPLE_NAME;

    public KeyMap(String keyFileName, String toolName, String BlankSampleName, int OUT_Q_CAPACITY) {
        TOOL_NAME = toolName;
        BLANK_SAMPLE_NAME = BlankSampleName;
        this.OUT_Q_CAPACITY = OUT_Q_CAPACITY;
        populateMap(keyFileName);
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
        sampleToCountMap = new ConcurrentHashMap<>(1000);
        BufferedReader content = null;
        try {
            if (inputFile.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile));
                content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"));
            } else {
                content = new BufferedReader(new FileReader(new File(inputFile)));
            }
            String line;
            while ((line = content.readLine()) != null && !line.isEmpty()) {
                String[] toks = line.split("\t");
                String flowcell = toks[0];
                String barcode = toks[2];
                String sample = toks[3];
                if(sample.compareTo(BLANK_SAMPLE_NAME) == 0) {
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
                        sampleToQueueMap.put(sample, new ArrayBlockingQueue<SampleBuffer>(OUT_Q_CAPACITY));
//                        sampleToCountMap.put(sample, 0L);
                    }
                } else {
                    ConcurrentHashMap<String, String> barcodeToSample = new ConcurrentHashMap<>();
                    barcodeToSample.put(barcode, sample);
                    keyMap.put(flowcell, barcodeToSample);
                    sampleToQueueMap.put(sample, new ArrayBlockingQueue<SampleBuffer>(OUT_Q_CAPACITY));
//                    sampleToCountMap.put(sample, 0L);
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

    public ConcurrentHashMap<String, BlockingQueue<SampleBuffer>> getSampleToQueueMap() {
        return sampleToQueueMap;
    }

    public ConcurrentHashMap<String, Long> getSampleToCountMap() {
        return sampleToCountMap;
    }

    public synchronized void addToSampleCount(String sample, long count) {
        Long stored = sampleToCountMap.get(sample); 
        if(stored == null) {
            stored = 0L;
        }
        sampleToCountMap.put(sample, stored+count);
    }
    
}
