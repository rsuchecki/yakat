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
package pileup2snps;

import java.util.ArrayList;
import java.util.regex.Pattern;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PileupPosition {

    private String referenceId;
    private int referencePosition;
    private char referenceBase;
    private ArrayList<SampleBases> sampleBases;

    public PileupPosition(String pileupLine, boolean IGNORE_SPLICING) {
        processPileup(pileupLine, IGNORE_SPLICING);
    }

    private void processPileup(String pileupLine, boolean IGNORE_SPLICING) {
        Pattern spliPattern = Pattern.compile("\t");
        String[] toks = spliPattern.split(pileupLine);
        sampleBases = new ArrayList<>((toks.length-3)/3);
        referenceId = toks[0];
        referencePosition = Integer.parseInt(toks[1]);
        referenceBase = toks[2].charAt(0);
        for (int i = 3; i < toks.length; i+=3) {
            sampleBases.add(new SampleBases(toks[i+1], referenceId, Integer.parseInt(toks[i]), IGNORE_SPLICING));
        }    
    }

    public String getReferenceId() {
        return referenceId;
    }

    public int getReferencePosition() {
        return referencePosition;
    }

    public char getReferenceBase() {
        return referenceBase;
    }
    
    public int getNumberOfSamples() {
        return sampleBases.size();
    }
    

}
