/*
 * Copyright 2016 rad.
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
package kmerextender;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class SeedSequences {

    private ArrayList<SeedSequence> seedSequences;
    private ArrayList<String> seedSequenceStrings;

    public SeedSequences(String fastaFileName, boolean expectingPurgedInput) {
        listOfSequencesFromFasta(fastaFileName, expectingPurgedInput);
    }

    public ArrayList<SeedSequence> getSeedSequences() {
        return seedSequences;
    }

    public ArrayList<String> getSeedSequenceStrings() {
        return seedSequenceStrings;
    }

    private void listOfSequencesFromFasta(String fastaFileName, boolean expectingPurgedInput) {
        seedSequences = new ArrayList<>();
        seedSequenceStrings = new ArrayList<>();
        File newFile = new File(fastaFileName);
        BufferedReader myData = null;
        try {
            String inputLine;
            myData = new BufferedReader(new FileReader(newFile));
            String id = "";
            StringBuilder seqBuilder = new StringBuilder();
            while ((inputLine = myData.readLine()) != null) {
                String line = inputLine.trim();
                if (line.startsWith(">")) {
                    if (seqBuilder.length() > 0) {
                        String seqString = seqBuilder.toString();
                        seedSequences.add(new SeedSequence(id, seqString));
                        if (!expectingPurgedInput) {
                            seedSequenceStrings.add(seqString);
                        }
                        seqBuilder = new StringBuilder();
                    }
                    id = line.substring(1); //get rid of ">"
                } else {
                    seqBuilder.append(line);
                }
            }
            if (id.isEmpty()) {
                System.err.println("Error reading FASTA [" + fastaFileName + "]. No identifier lines? Terminating... ");
                System.exit(1);
            } else {
//                String identifierString[] = id.split(" ");
                String seqString = seqBuilder.toString();
                seedSequences.add(new SeedSequence(id, seqString));
                if (!expectingPurgedInput) {
                    seedSequenceStrings.add(seqString);
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + newFile.getName(), this.getClass().getName());
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (myData != null) {
                    myData.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
}
