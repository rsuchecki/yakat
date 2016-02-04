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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class InputReaderProducer implements Runnable {

    private BlockingQueue queue; //Either
    private HashMap<Integer, BlockingQueue> kSizeToQueue; //or

    private ArrayList<Integer> kValues;

    private Integer KMER_LENGTH; // ignored if <0 but not if null
    private ArrayList<String> inputFiles;
    private GuessedInputFormat guessedInputFormat;
    private final int READER_BUFFER_SIZE = 8192;
//    private MerMap map; //almost dummy object used to pass out of mem error up the chain
    //PUT READ LINES ON QUEQE AS SOON AS ONE OF THE FOLLOWING IS REACHED
//    private final int FASTA_RECORD_LEN = 8192;    
//    private final int FASTA_BUFFER_SIZE = 128; //THAT MANY FASTA RECORDS
    private int FASTQ_BUFFER_SIZE = 1024; //THAT MANY FASTQ RECORDS
    private int KMER_BUFFER_SIZE = 8192; // //THAT MANY KMERS 
//    private final int KMER_REPORTING_MULTIPLY = 2; //nice to use 2 if KMER_BUFFER_SIZE is a power of2 or 10 if it is a power of 10 
    private final String TOOL_NAME;
//    private String RECORD_NAME = "kmers";

    public enum GuessedInputFormat {
        KMERS, FASTA_SE_ONE_LINE, FASTA_PE_ONE_LINE, FASTQ_SE_ONE_LINE, FASTQ_PE_ONE_LINE, FASTA, FASTQ, UNSUPPORTED_OR_UNRECOGNIZED;
    }

    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, Integer k, String toolName, String RECORD_NAME, int RECORD_BUFFER_SIZE) {
        this.queue = queue;
        this.inputFiles = inputFiles;
        KMER_LENGTH = k;
        TOOL_NAME = toolName;
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
//        this.RECORD_NAME = RECORD_NAME;
    }

    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, Integer k, MerMap map, String toolName) {
        this.queue = queue;
        this.inputFiles = inputFiles;
        KMER_LENGTH = k;
//        if (map != null) {
//            this.map = map;
//        }
        TOOL_NAME = toolName;
    }

    /**
     * Used for kmer extender
     *
     * @param kSizeToQueue
     * @param kValues
     * @param inputFiles
     * @param toolName
     */
    public InputReaderProducer(HashMap<Integer, BlockingQueue> kSizeToQueue, ArrayList<Integer> kValues,
        ArrayList<String> inputFiles, int RECORD_BUFFER_SIZE, String toolName) {
        if (kSizeToQueue.size() == 1) {
            this.queue = kSizeToQueue.values().iterator().next();
            this.KMER_LENGTH = kSizeToQueue.keySet().iterator().next();
        } else {
            this.kSizeToQueue = kSizeToQueue;
            this.kValues = kValues;
        }
        this.inputFiles = inputFiles;

//        if (map != null) {
//            this.map = map;
//        }
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        TOOL_NAME = toolName;
    }

    /**
     * Execute InputReaderProducer thread
     */
    @Override
    public void run() throws OutOfMemoryError {
        int numTestLines = 8;
        BufferedReader content = null;
        try {
            if (inputFiles == null || inputFiles.isEmpty()) { //READSTDIN
                inputFiles = new ArrayList<>();
                inputFiles.add("-");
            }
            for (String inputFile : inputFiles) {
                if (inputFile.equals("-")) { //READSTDIN
                    content = new BufferedReader(new InputStreamReader(System.in), READER_BUFFER_SIZE);
                } else if (inputFile.endsWith(".gz")) {// reading kmers from a compressed file
                    InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile), READER_BUFFER_SIZE);
                    content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
                } else {
                    content = new BufferedReader(new FileReader(new File(inputFile)), READER_BUFFER_SIZE);//reading kmers from a text file
                }
                String line1;
                //TAKE UP TO 4 FIRST LINES AND TRY TO GUESS THE INPUT FORMAT
                ArrayList<String> testLines = new ArrayList<>(numTestLines);
                for (int i = 0; i < numTestLines; i++) {
                    if ((line1 = content.readLine()) != null && !line1.isEmpty()) {
                        if (!line1.trim().isEmpty()) {
                            testLines.add(line1.trim());
                        }
                    }
                }
                //IF FORMAT SUPPORT NOT IMPLEMENTED OR UNRECOGNZED 
                guessedInputFormat = guessInputFormat(testLines);
                Reporter.report("[INFO]", "Input format guessed: " + guessedInputFormat.toString(), TOOL_NAME);
                if (guessedInputFormat == GuessedInputFormat.UNSUPPORTED_OR_UNRECOGNIZED) {
                    Reporter.report("[ERROR]", "Unrecognized or unsupported input, terminating...", TOOL_NAME);
                    putOnQueue(new ArrayList<String>()); //OTHERWISE CONSUMER THREAD WILL KEEP GOING
                    System.exit(1);
                }

                String line;
                if (guessedInputFormat == GuessedInputFormat.KMERS || (guessedInputFormat == GuessedInputFormat.FASTQ_PE_ONE_LINE && !kMerIsSet()) || (guessedInputFormat == GuessedInputFormat.FASTQ_SE_ONE_LINE && !kMerIsSet())) {
                    //READ KMERS
                    readKmers(content, testLines);
                } else if (!kMerIsSet()) {
                    Reporter.report("[ERROR]", "Fatal error, k-mer lenght must be specified for input other than a list of k-mers, terminating!", TOOL_NAME);
                    putOnQueue(new ArrayList<String>(0)); //OTHERWISE CONSUMER THREAD WILL KEEP GOING
                    System.exit(1);
                } else if (guessedInputFormat == GuessedInputFormat.FASTA) {
                    StringBuilder sb = new StringBuilder();
                    //Don't forget the test lines 
                    for (String string : testLines) {
                        sb = processFastaLine(string, sb);
                    }
                    //read subsequent lines
                    while ((line = content.readLine()) != null && !line.isEmpty()) {
                        sb = processFastaLine(line, sb);
                    }
                    //make sure last record is put on the queue
                    processFastaLine(">", sb);
                } else if (guessedInputFormat == GuessedInputFormat.FASTQ) {
                    //READ FASTQ
                    readFastq(content, testLines);
                }
            }
//            queue.put(new ArrayList<String>()); //TELLS CONSUMERS, NO MORE DATA
            putOnQueue(new ArrayList<String>(0));
        } catch (OutOfMemoryError err) {
//            if (map == null) {
            Reporter.report("[ERROR]", "Out of memory error!", TOOL_NAME);
//            } else {
//                map.setOutOfMemory();
//            }
            try {
//                queue.put(new ArrayList<String>());
                putOnQueue(new ArrayList<String>(0));
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
        } catch (InterruptedException ex) {
            System.err.println(ex.getMessage());
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

    private void readKmers(BufferedReader content, ArrayList<String> testLines) 
        throws InterruptedException, IOException {
        ArrayList<String> bufferList = new ArrayList<>(KMER_BUFFER_SIZE);
        bufferList.addAll(testLines);
//                    queue.put(testLines);
        long kmerCount = 0L;
        long reportThreshold = (long) KMER_BUFFER_SIZE;
        String line;
        while ((line = content.readLine()) != null && !line.isEmpty()) {
            if (bufferList.size() == KMER_BUFFER_SIZE) {
//                            queue.put(bufferList);
                putOnQueue(bufferList);
                kmerCount += KMER_BUFFER_SIZE;
                bufferList = new ArrayList<>();
                if (kmerCount % reportThreshold == 0) {
                    if (reportThreshold < 1e9) {
//                                    reportThreshold *= 10;
//                                    reportThreshold *= KMER_REPORTING_MULTIPLY;
                        reportThreshold <<= 1; // *= 2
                    }
                    Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInputFormat.toString() + " read-in so far", TOOL_NAME);
                }
            }
            bufferList.add(line);
//                        queue.put(line);
        }
        kmerCount += bufferList.size();
        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInputFormat.toString() + " read-in", TOOL_NAME);
//                    queue.put(bufferList);
        putOnQueue(bufferList);
    }

    private void readFastq(BufferedReader content, ArrayList<String> testLines) 
        throws InterruptedException, IOException {
        long fastqCount = 0L;
        long reportThreshold = (long) FASTQ_BUFFER_SIZE;
        int countFastqLine = 1; //WE ONLY WANT THE NUCL SEQ FROM FASTQ NOT THE OTHER THREE LINES OF EACH RECORD
        ArrayList<String> bufferList = new ArrayList<>(FASTQ_BUFFER_SIZE);
        for (int i = 1; i < testLines.size(); i += 4) {
            bufferList.add(testLines.get(i));
        }
        String line;
        while ((line = content.readLine()) != null && !line.isEmpty()) {
            if (countFastqLine++ == 2) {
                if (bufferList.size() == FASTQ_BUFFER_SIZE) {
//                                queue.put(bufferList);
                    putOnQueue(bufferList);
                    bufferList = new ArrayList<>();
                    fastqCount += FASTQ_BUFFER_SIZE;
                    if (fastqCount % reportThreshold == 0) {
                        if (reportThreshold < 1e9) {
//                                    reportThreshold *= 10;
//                                    reportThreshold *= KMER_REPORTING_MULTIPLY;
                            reportThreshold <<= 1; // *= 2
                        }
                        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(fastqCount) + " " + guessedInputFormat.toString() + " read-in so far", TOOL_NAME);
                    }
                }
                bufferList.add(line);
            }
            if (countFastqLine % 4 == 0) {
                countFastqLine = 0;
            }
        }
//                    queue.put(bufferList);
        putOnQueue(bufferList);
        fastqCount += bufferList.size();
        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(fastqCount) + " " + guessedInputFormat.toString() + " read-in", TOOL_NAME);
    }

    private StringBuilder processFastaLine(String line, StringBuilder sb) throws InterruptedException {
        if (line.startsWith(">")) {
            if (sb.length() > 0) {
//                bufferList.add(sb.toString());
//                if(bufferList.size() >= FASTA_BUFFER_SIZE || sb.length() >= FASTA_RECORD_LEN || line.endsWith(">")) {
//                    queue.put(bufferList);
//                    bufferList = new ArrayList<>(FASTA_BUFFER_SIZE);
//                }

                ArrayList<String> wrapper = new ArrayList<>(1);
                wrapper.add(sb.toString());
//                queue.put(wrapper);
                putOnQueue(wrapper);
            }
            return new StringBuilder();
        } else {
            return sb.append(line);
        }
    }

    /**
     * Attempt to guess the input format and set KMER_LNGTH if KMERS input
     *
     * @param testLines
     * @return
     */
    private GuessedInputFormat guessInputFormat(ArrayList<String> testLines) {
        if (testLines.isEmpty()) {
            return GuessedInputFormat.UNSUPPORTED_OR_UNRECOGNIZED;
        }
        if (!testLines.get(0).startsWith("@") && !testLines.get(0).startsWith(">")) {
            String[] split0 = testLines.get(0).split("\t| ");
            String[] split1 = testLines.get(1).split("\t| ");
            if (split0[0].matches("^[A|T|C|G]+$") && split1[0].matches("^[A|T|C|G]+$") && split0[0].length() == split1[0].length()) {
                setKmerLength(split0[0].trim().length());
                return GuessedInputFormat.KMERS;
            } else {
                return GuessedInputFormat.UNSUPPORTED_OR_UNRECOGNIZED;
            }
        }

        if (testLines.size() >= 4) {
            if (testLines.get(0).startsWith("@") && testLines.get(2).startsWith("+")) {
                return GuessedInputFormat.FASTQ;
            }
        }
        if (testLines.size() > 1) {
            if (testLines.get(0).startsWith(">") && !testLines.get(1).startsWith(">")) {
                return GuessedInputFormat.FASTA;
            }
        }

        String[] split = testLines.get(0).split("\t");
        if (split.length == 2 && split[0].startsWith(">")) {
            return GuessedInputFormat.FASTA_SE_ONE_LINE;
        } else if (split.length == 4) {
            if (split[0].startsWith(">") && split[2].startsWith(">")) {
                return GuessedInputFormat.FASTA_PE_ONE_LINE;
            }
            if (split[0].startsWith("@") && split[2].startsWith("+")) {
                return GuessedInputFormat.FASTQ_SE_ONE_LINE;
            }
        } else if (split.length == 8) {
            if (split[0].startsWith("@") && split[2].startsWith("+") && split[4].startsWith("@") && split[6].startsWith("+")) {
                return GuessedInputFormat.FASTQ_PE_ONE_LINE;
            }
        }

        return GuessedInputFormat.UNSUPPORTED_OR_UNRECOGNIZED;
    }

    public GuessedInputFormat getGuessedInputFormat() {
        return guessedInputFormat;
    }

    private void setKmerLength(int k) {
        this.KMER_LENGTH = k;
    }

    public Integer getKmerLength() {
        return KMER_LENGTH;
    }

    private boolean kMerIsSet() {
        return (KMER_LENGTH != null && KMER_LENGTH > 0) || (kValues != null && !kValues.isEmpty());
    }

    private void putOnQueue(ArrayList<String> list) throws InterruptedException {
        if (queue != null) {
            queue.put(list);
        } else {
            for (Integer k : kValues) {
                kSizeToQueue.get(k).put(list);
            }
        }
    }
}
