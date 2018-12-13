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

import java.util.ArrayList;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class PerSampleBuffer {
//    private String sampleId;
    private final String sampleId;
    private final String outFileName;
    private boolean appendToOutputFile;
    private ArrayList<String> bufferedRecords;

    public PerSampleBuffer(String sampleId, String outFileName, int bufferSize) {
        this.sampleId = sampleId;
        this.outFileName = outFileName;
        this.bufferedRecords = new ArrayList<>(bufferSize);
    }

    public PerSampleBuffer() {
        this.sampleId = null;
        this.outFileName = null;
    }
    
    

//    public PerSampleBuffer resetBuffer(int bufferSize) {
//        this.bufferedRecords = new ArrayList<>(bufferSize);        
//        return this;
//    }
    
    public String getSampleId() {
        return sampleId;
    }
//
//    public String getOutFileName() {
//        return outFileName;
//    }
//
//    public boolean isAppendToOutputFile() {
//        return appendToOutputFile;
//    }
//
//    public ArrayList<String> getBufferedRecords() {
//        return bufferedRecords;
//    }

    public int size() {
        return bufferedRecords.size();
    }
    
    public void add(String record) {
        bufferedRecords.add(record);
    }
    
    public boolean isEmpty() {
        return bufferedRecords == null || bufferedRecords.isEmpty();
    }
     
    public String remove(int index) {
        return bufferedRecords.remove(index);
    }
    
    
}
