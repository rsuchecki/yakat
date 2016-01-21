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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class FastaReader {

    /**
     * Reads in a fasta file and generates a HashMap of Sequence objects input
     * fasta sequences can be broken over multiple lines
     *
     * @param fastaFileName
     * @return
     */
    public static ArrayList<Sequence> listOfSequencesFromFasta(String fastaFileName, Set<String> requiredContigIds) {
        ArrayList<Sequence> sequencesList = new ArrayList<>();
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
                        String identifierString[] = id.split(" ");
                        String key = identifierString[0];
                        if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                            sequencesList.add(new Sequence(id, seqBuilder.toString()));
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
                String identifierString[] = id.split(" ");
                String key = identifierString[0];
                if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                    sequencesList.add(new Sequence(id, seqBuilder.toString()));
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + newFile.getName(), getClassName());
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
        return sequencesList;
    }
    

    /**
     * Reads in a fasta file and generates a HashMap of Sequence objects input
     * fasta sequences can be broken over multiple lines
     *
     * @param fastaFileName
     * @return
     */
    public static HashMap<String, Sequence> hashMapOfSequencesFromFasta(String fastaFileName, Set<String> requiredContigIds) {
        HashMap<String, Sequence> sequencesMap = new HashMap<>();
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
                        String identifierString[] = id.split(" ");
                        String key = identifierString[0];
                        if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                            Sequence previous = sequencesMap.put(key, new Sequence(id, seqBuilder.toString()));
                            if (previous != null) {
                                System.err.println("Fatal error: non unique sequence identifier: " + key + " in " + fastaFileName);
                                System.exit(1);
                            }
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
                String identifierString[] = id.split(" ");
                String key = identifierString[0];
                if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                    Sequence previous = sequencesMap.put(key, new Sequence(id, seqBuilder.toString()));
                    if (previous != null) {
                        System.err.println("Fatal error: non unique sequence identifier: " + key + " in " + fastaFileName);
                        System.exit(1);
                    }
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + newFile.getName(), getClassName());
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
        return sequencesMap;
    }

    public static HashMap<String, Integer> hashMapOfSequenceLengthsFromFasta(String fastaFileName, Set<String> requiredContigIds) {
        HashMap<String, Integer> sequencesMap = new HashMap<>();
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
                        String identifierString[] = id.split(" ");
                        String key = identifierString[0];
                        if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                            Integer previous = sequencesMap.put(key, seqBuilder.length());
                            if (previous != null) {
                                System.err.println("Fatal error: non unique sequence identifier: " + key + " in " + fastaFileName);
                                System.exit(1);
                            }
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
                String identifierString[] = id.split(" ");
                String key = identifierString[0];
                if (requiredContigIds == null || requiredContigIds.contains(key)) { //addedd to preserve memory only!
                    Integer previous = sequencesMap.put(key, seqBuilder.length());
                    if (previous != null) {
                        System.err.println("Fatal error: non unique sequence identifier: " + key + " in " + fastaFileName);
                        System.exit(1);
                    }
                }
            }

        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + newFile.getName(), getClassName());
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
        return sequencesMap;
    }

    private static String getClassName() {
        return MethodHandles.lookup().lookupClass().toString().replaceFirst("class ", "");
    }
}
