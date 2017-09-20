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
import java.util.concurrent.BlockingQueue;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class InputReaderProducer implements Runnable {

    private final BlockingQueue queue; //Either
//    private HashMap<Integer, BlockingQueue> kSizeToQueue; //or

    private ArrayList<Integer> kLengths;

//    private Integer KMER_LENGTH; // ignored if <0 but not if null
    private ArrayList<String> inputFiles;
    private InFormat guessedInFormat;
    private boolean ASSUME_GENERIC_ONE_RECORD_PER_LINE_FORMAT;
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
    private boolean useLabelledBuffers; //so far only used for snpmers module
    private int REPORTING_SHIFT = 0;

    public enum InFormat {
        KMERS, PER_SAMPLE_KMERS, FASTA_SE_ONE_LINE, FASTA_PE_ONE_LINE, FASTQ_SE_ONE_LINE, FASTQ_PE_ONE_LINE, FASTQ_PE_WITH_INDEX_ONE_LINE, FASTA, FASTQ,
        MPILEUP, ONE_RECORD_PER_LINE, UNSUPPORTED_OR_UNRECOGNIZED, EMPTY;
    }

    /**
     * Used for reading pileup
     *
     * @param queue
     * @param inputFile
     * @param k
     * @param toolName
     * @param RECORD_BUFFER_SIZE
     */
    public InputReaderProducer(BlockingQueue queue, String inputFile, Integer k, String toolName, int RECORD_BUFFER_SIZE) {
        this.queue = queue;
        if (inputFile != null) {
            this.inputFiles = new ArrayList<>(1);
            inputFiles.add(inputFile);
        }
//        KMER_LENGTH = k;
        kLengths = new ArrayList<>(1);
        kLengths.add(k);
        TOOL_NAME = toolName;
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
//        this.RECORD_NAME = RECORD_NAME;
    }

    /**
     * Used by splitgbs
     *
     * @param queue
     * @param inputFiles
     * @param toolName
     * @param RECORD_NAME
     * @param RECORD_BUFFER_SIZE
     */
    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, String toolName, String RECORD_NAME, int RECORD_BUFFER_SIZE) {
        this.queue = queue;
        this.inputFiles = inputFiles;
        TOOL_NAME = toolName;
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
//        this.RECORD_NAME = RECORD_NAME;
    }

    /**
     * Used by idmatch
     *
     * @param queue
     * @param inputFile
     * @param toolName
     * @param RECORD_NAME
     * @param RECORD_BUFFER_SIZE
     */
    public InputReaderProducer(BlockingQueue queue, String inputFile, String toolName, String RECORD_NAME, int RECORD_BUFFER_SIZE, boolean assumeInputFormat) {
        this.queue = queue;
        this.inputFiles = new ArrayList<>();
        inputFiles.add(inputFile);
        TOOL_NAME = toolName;
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        ASSUME_GENERIC_ONE_RECORD_PER_LINE_FORMAT = assumeInputFormat;
//        this.RECORD_NAME = RECORD_NAME;
    }

    /**
     *
     * @param queue
     * @param inputFiles
     * @param k
     * @param toolName
     * @param RECORD_NAME
     * @param RECORD_BUFFER_SIZE
     */
    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, Integer k, String toolName, String RECORD_NAME, int RECORD_BUFFER_SIZE) {
        this.queue = queue;
        this.inputFiles = inputFiles;
//        KMER_LENGTH = k;
        kLengths = new ArrayList<>(1);
        kLengths.add(k);
        TOOL_NAME = toolName;
        KMER_BUFFER_SIZE = RECORD_BUFFER_SIZE;
        FASTQ_BUFFER_SIZE = RECORD_BUFFER_SIZE;
//        this.RECORD_NAME = RECORD_NAME;
    }

    /**
     * used by (unfinished) kmer-match
     *
     * @param queue
     * @param inputFiles
     * @param k
     * @param map
     * @param toolName
     */
    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, Integer k, MerMap map, String toolName) {
        this.queue = queue;
        this.inputFiles = inputFiles;
        kLengths = new ArrayList<>(1);
        kLengths.add(k);
//        KMER_LENGTH = k;
//        if (map != null) {
//            this.map = map;
//        }
        TOOL_NAME = toolName;
    }

    /**
     * Used by SnpMers
     *
     * @param queue
     * @param inputFiles
     * @param useLabelledBuffers
     * @param k
     * @param map
     * @param toolName
     */
    public InputReaderProducer(BlockingQueue queue, ArrayList<String> inputFiles, boolean useLabelledBuffers, int reportingShift,
            String toolName, int KMER_BUFFER_SIZE) {
        this.queue = queue;
        this.inputFiles = inputFiles;
        this.useLabelledBuffers = useLabelledBuffers;
        this.TOOL_NAME = toolName;
        this.KMER_BUFFER_SIZE = KMER_BUFFER_SIZE;
        this.REPORTING_SHIFT = reportingShift;
    }

    /**
     * Used for kmer extender
     *
     * @param queue
     * @param kSizes
     * @param RECORD_BUFFER_SIZE
     * @param inputFiles
     * @param toolName
     */
//    public InputReaderProducer(HashMap<Integer, BlockingQueue> kSizeToQueue, ArrayList<Integer> kValues,
    public InputReaderProducer(BlockingQueue queue, ArrayList<Integer> kSizes,
            ArrayList<String> inputFiles, int RECORD_BUFFER_SIZE, String toolName) {
//        if (kSizeToQueue.size() == 1) {
        this.queue = queue; // kSizeToQueue.values().iterator().next();
//            this.KMER_LENGTH = kSizeToQueue.keySet().iterator().next();
//        } else {
//            this.kSizeToQueue = kSizeToQueue;
//        }
        this.kLengths = kSizes;
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
            Object empty = new ArrayList<>(0);
            for (String inputFile : inputFiles) {
                if (inputFile.equals("-")) { //READSTDIN
                    content = new BufferedReader(new InputStreamReader(System.in, "UTF-8"), READER_BUFFER_SIZE);
                } else if (inputFile.endsWith(".gz")) {// reading kmers from a compressed file
                    InputStream gzipStream = new GZIPInputStream(new FileInputStream(inputFile), READER_BUFFER_SIZE);
                    content = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
                } else {
                    content = new BufferedReader(new FileReader(new File(inputFile)), READER_BUFFER_SIZE);//reading from a text file
                }

                //Non-UNIX line endings detected
                content.mark(65535);
                char[] header = new char[65535];
                content.read(header);
                for (int i = 0; i < header.length; i++) {
                    char c = header[i];
                    if (c == '\r') {
                        Reporter.report("[FATAL]", "Non-UNIX line ending char detected (\\r), terminating...", TOOL_NAME);
                        putOnQueue(new ArrayList<>()); //OTHERWISE CONSUMER THREADS WILL KEEP GOING
                        System.exit(1);
                    }
                }
                content.reset();

                String line1;
                //TAKE UP TO 8 FIRST LINES AND TRY TO GUESS THE INPUT FORMAT
                ArrayList<String> testLines = new ArrayList<>(numTestLines);
                for (int i = 0; i < numTestLines; i++) {
                    if ((line1 = content.readLine()) != null && !line1.isEmpty()) {
                        if (!line1.trim().isEmpty()) {
                            testLines.add(line1.trim());
                        }
                    }
                }
                //IF FORMAT SUPPORT NOT IMPLEMENTED OR UNRECOGNZED 
                if (ASSUME_GENERIC_ONE_RECORD_PER_LINE_FORMAT) {
                    guessedInFormat = InFormat.ONE_RECORD_PER_LINE;
                    Reporter.report("[INFO]", "Input format forced: " + guessedInFormat.toString(), TOOL_NAME);
                } else {
                    guessedInFormat = guessInputFormat(testLines);
                    Reporter.report("[INFO]", "Input format guessed: " + guessedInFormat.toString(), TOOL_NAME);
                }
                if (guessedInFormat == InFormat.EMPTY) {
                    Reporter.report("[WARNING]", "Empty input file/stream...", TOOL_NAME);
                    putOnQueue(new ArrayList<>()); //OTHERWISE CONSUMER THREADS WILL KEEP GOING                    
                }
                if (guessedInFormat == InFormat.UNSUPPORTED_OR_UNRECOGNIZED) {
                    Reporter.report("[FATAL]", "Unrecognized or unsupported input, terminating...", TOOL_NAME);
                    putOnQueue(new ArrayList<>()); //OTHERWISE CONSUMER THREADS WILL KEEP GOING
                    System.exit(1);
                }

                String line;
                if (guessedInFormat == InFormat.PER_SAMPLE_KMERS || (guessedInFormat == InFormat.KMERS && useLabelledBuffers)) {
                    //Package k-mers so that the consumer will know what sample/dataset they came from based on the sample label                    
                    readKmersPerSample(content, testLines, inputFile);
                    empty = new LabelledInputBuffer("", new ArrayList()); //making sure that the type of "done reading" flag matches the input type
                } else if (guessedInFormat == InFormat.KMERS
                        || (guessedInFormat == InFormat.FASTQ_PE_ONE_LINE && !kMerIsSet())
                        || (guessedInFormat == InFormat.FASTQ_PE_WITH_INDEX_ONE_LINE && !kMerIsSet())
                        || (guessedInFormat == InFormat.FASTQ_SE_ONE_LINE && !kMerIsSet())
                        || (guessedInFormat == InFormat.MPILEUP)
                        || (guessedInFormat == InFormat.ONE_RECORD_PER_LINE)) {
                    //READ KMERS (or any other input that can simply be processed line by line)
                    readLines(content, testLines);
                } else if (!kMerIsSet()) {
                    Reporter.report("[FATAL]", "Fatal error, k-mer length must be specified for input other than a list of k-mers, terminating!", TOOL_NAME);
                    putOnQueue(new ArrayList<String>(0)); //OTHERWISE CONSUMER THREAD WILL KEEP GOING
                    System.exit(1);
                } else if (guessedInFormat == InFormat.FASTA) {
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
                } else if (guessedInFormat == InFormat.FASTQ) {
                    //READ FASTQ
                    readSequenceFromFastq(content, testLines);
                }
            }
//            queue.put(new ArrayList<String>()); //TELLS CONSUMERS, NO MORE DATA
            putOnQueue(empty);
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

    private void readLines(BufferedReader content, ArrayList<String> testLines)
            throws InterruptedException, IOException {
        ArrayList<String> bufferList = new ArrayList<>(KMER_BUFFER_SIZE);
        bufferList.addAll(testLines);
        long kmerCount = 0L;
        long reportThreshold = (long) KMER_BUFFER_SIZE;
        reportThreshold <<= REPORTING_SHIFT;
        String line;
        while ((line = content.readLine()) != null && !line.isEmpty()) {
            if (bufferList.size() == KMER_BUFFER_SIZE) {
                putOnQueue(bufferList);
                kmerCount += KMER_BUFFER_SIZE;
                bufferList = new ArrayList<>();
                if (kmerCount % reportThreshold == 0) {
                    if (reportThreshold < 1e8) {
//                                    reportThreshold *= 10;
//                                    reportThreshold *= KMER_REPORTING_MULTIPLY;
                        reportThreshold <<= 1; // *= 2
                    }
                    Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInFormat.toString() + " read-in so far", TOOL_NAME);
                }
            }
            bufferList.add(line);
        }
        kmerCount += bufferList.size();
        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInFormat.toString() + " read-in", TOOL_NAME);
        putOnQueue(bufferList);
    }

    /**
     * Used when we need to know where given lines have come from. This can be
     * for example a sample label or a file name
     *
     * @param content
     * @param testLines
     * @throws InterruptedException
     * @throws IOException
     */
    private void readKmersPerSample(BufferedReader content, ArrayList<String> testLines, String inputFile)
            throws InterruptedException, IOException {
//        Reporter.report("[WARNING]", "readKmersPerSample() needs to be tested, particularily with multiple input files and or multiple samples per file", TOOL_NAME);
        String sampleLabel = null;
        if (!inputFile.equals("-")) { //IF NOT STDIN
            sampleLabel = inputFile;  //Assume filename is the label
        }
        ArrayList<String[]> bufferList = new ArrayList<>(KMER_BUFFER_SIZE);
        long kmerCount = 0L;
        long reportThreshold = (long) KMER_BUFFER_SIZE;
        reportThreshold <<= REPORTING_SHIFT;

        //PROCESS TESTLINES
        for (String line : testLines) {
            String[] toks = line.split("\t");
            if (toks.length == 1) {
                sampleLabel = line;
                if (!bufferList.isEmpty()) {
                    putOnQueue(new LabelledInputBuffer(sampleLabel, bufferList));
                    kmerCount += bufferList.size();
                    bufferList = new ArrayList<>();
                }
            } else {
                bufferList.add(toks);
            }
        }
        String line;
        //PROCESS REMINDER OF THE INPUT
        while ((line = content.readLine()) != null && !line.isEmpty()) {
            String[] toks = line.split("\t");
            if (toks.length == 1 || bufferList.size() == KMER_BUFFER_SIZE) {
                if (!bufferList.isEmpty()) {
                    putOnQueue(new LabelledInputBuffer(sampleLabel, bufferList));
                    kmerCount += bufferList.size();
                    bufferList = new ArrayList<>();
                    if (kmerCount % reportThreshold == 0) {
                        if (reportThreshold < 1e9) {
                            reportThreshold <<= 1; // *= 2
                        }
                        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInFormat.toString() + " read-in so far", TOOL_NAME);
                    }
                }
                if (toks.length == 1) {
                    Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInFormat.toString() + " read-in for sample " + sampleLabel, TOOL_NAME);
                    kmerCount = 0;
                    reportThreshold = (long) KMER_BUFFER_SIZE;
                    reportThreshold <<= REPORTING_SHIFT;
                    sampleLabel = toks[0];
                } else {
                    bufferList.add(toks);
                }
            } else {
                bufferList.add(toks);
            }
        }
        kmerCount += bufferList.size();
        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(kmerCount) + " " + guessedInFormat.toString() + " read-in for sample " + sampleLabel, TOOL_NAME);
//                    queue.put(bufferList);
        putOnQueue(new LabelledInputBuffer(sampleLabel, bufferList));
    }

    private void readSequenceFromFastq(BufferedReader content, ArrayList<String> testLines)
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
                        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(fastqCount) + " " + guessedInFormat.toString() + " read-in so far", TOOL_NAME);
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
        Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(fastqCount) + " " + guessedInFormat.toString() + " read-in", TOOL_NAME);
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
    private InFormat guessInputFormat(ArrayList<String> testLines) {
        if (testLines.isEmpty()) {
            return InFormat.EMPTY;
        }
        if (!testLines.get(0).startsWith("@") && !testLines.get(0).startsWith(">") & testLines.size() > 1) {
            String[] split0 = testLines.get(0).split("\t| ");
            String[] split1 = testLines.get(1).split("\t| ");
            if (split0[0].matches("^[A|T|C|G]+$") && split1[0].matches("^[A|T|C|G]+$") && split0[0].length() == split1[0].length()) {
//                setKmerLength(split0[0].trim().length());
                addKmerLength(split0[0].trim().length());
                return InFormat.KMERS;
//            } else {
//                return InFormat.UNSUPPORTED_OR_UNRECOGNIZED;
            }
            if (testLines.size() > 2) {
                String[] split2 = testLines.get(2).split("\t| ");
                if (split0.length == 1 && split1[0].matches("^[A|T|C|G]+$") && split2[0].matches("^[A|T|C|G]+$") && split1[0].length() == split2[0].length()) {
                    addKmerLength(split1[0].trim().length());
                    return InFormat.PER_SAMPLE_KMERS;
                }
            }
        }

        if (testLines.size() >= 4) {
            if (testLines.get(0).startsWith("@") && (testLines.get(2).startsWith("+") || testLines.get(2).equals(testLines.get(0)))) {
                return InFormat.FASTQ;
            }
        }
        if (testLines.size() > 1) {
            if (testLines.get(0).startsWith(">") && !testLines.get(1).startsWith(">")) {
                return InFormat.FASTA;
            }
        }

        String[] split = testLines.get(0).split("\t");
        if (split.length == 2 && split[0].startsWith(">")) {
            return InFormat.FASTA_SE_ONE_LINE;
        } else if (split.length == 4) {
            if (split[0].startsWith(">") && split[2].startsWith(">")) {
                return InFormat.FASTA_PE_ONE_LINE;
            }
            if (split[0].startsWith("@") && split[2].startsWith("+")) {
                return InFormat.FASTQ_SE_ONE_LINE;
            }
        } else if (split.length == 8) {
            if (split[0].startsWith("@") && (split[2].startsWith("+") || split[2].startsWith("@"))
                    && split[4].startsWith("@") && (split[6].startsWith("+") || split[6].startsWith("@"))) {
                return InFormat.FASTQ_PE_ONE_LINE;
            }
        } else if (split.length == 12) {
            if (split[0].startsWith("@") && (split[2].startsWith("+") || split[2].startsWith("@"))
                    && split[4].startsWith("@") && (split[6].equals("+") || split[6].startsWith("@"))
                    && split[8].startsWith("@") && (split[10].equals("+") || split[10].startsWith("@"))) {
                return InFormat.FASTQ_PE_WITH_INDEX_ONE_LINE;
            }
        }

        //MPILEUP ref refPos refBase then 3 cols per sample
        if (split.length > 5 && ((split.length - 3) % 3 == 0) && split[2].length() == 1) {
            try {
                Integer.parseInt(split[1]);
                boolean mpileup = true;
                //looks like mpileup so far
                for (int i = 3; i < split.length; i += 3) {
                    try {
                        Integer.parseInt(split[i]);
                    } catch (NumberFormatException e) {
                        mpileup = false;
                    }
                }
                if (mpileup) {
                    return InFormat.MPILEUP;
                }
            } catch (NumberFormatException e) {
            }
        }

        return InFormat.UNSUPPORTED_OR_UNRECOGNIZED;
    }

    public InFormat getGuessedInputFormat() {
        return guessedInFormat;
    }

//    private void setKmerLength(int k) {
//        this.KMER_LENGTH = k;
//    }
    private synchronized void addKmerLength(int k) {
        if (kLengths == null) {
            kLengths = new ArrayList<>();
        }
        //if list contains just a dummy entrance (0) - replace it
        if (kLengths.size() == 1 && kLengths.get(0) == 0) {
            kLengths.remove(0);
        }
//        System.err.print("\n"+kLengths.size() + " kSizes:");
//            for (Integer kSize : kLengths) {
//                System.err.print(" " + kSize);
//            }
        if (!kLengths.contains(k)) {
            kLengths.add(k);
            Reporter.report("[INFO]", "Inreader adding k=" + k, TOOL_NAME);
        }

    }

    public ArrayList<Integer> getKmerLengths() {
        return kLengths;
    }
//    public Integer getKmerLength() {
//        return KMER_LENGTH;
//    }

    private boolean kMerIsSet() {
//        return (KMER_LENGTH != null && KMER_LENGTH > 0) || (kValues != null && !kValues.isEmpty());
        return kLengths != null && !kLengths.isEmpty() && kLengths.get(0) > 0;
    }

//    private void putOnQueue(ArrayList<String> list) throws InterruptedException {
    private void putOnQueue(Object list) throws InterruptedException {
//        if (queue != null) {

        queue.put(list);
//        } else {
//            for (Integer k : kLengths) {
//                kSizeToQueue.get(k).put(list);
//            }
//        }
    }
}
