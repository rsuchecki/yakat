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
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SampleBuffer {
//    private String sampleId;
    private final String sampleId;
    private final String outFileName;
    private boolean appendToOutputFile;
    private ArrayList<String> bufferedRecords;

    public SampleBuffer(String sampleId, String outFileName, ArrayList<String> bufferedRecords) {
        this.sampleId = sampleId;
        this.outFileName = outFileName;
        this.bufferedRecords = bufferedRecords;
    }

    public String getSampleId() {
        return sampleId;
    }

    public String getOutFileName() {
        return outFileName;
    }

    public boolean isAppendToOutputFile() {
        return appendToOutputFile;
    }

    public ArrayList<String> getBufferedRecords() {
        return bufferedRecords;
    }

    
    
    
    
}
