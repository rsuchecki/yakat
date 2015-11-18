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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.zip.GZIPOutputStream;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class FileWriterConsumer implements Runnable {

    private final BlockingQueue<ArrayList<String>> outputQueue;
    private final int BUFFER_SIZE = 8192; // //THAT MANY KMERS 
    private String RECORD_NAME = "lines";
    private final String TOOL_NAME;
    private final String FILE_NAME;

    public FileWriterConsumer(String fileName, BlockingQueue<ArrayList<String>> outputQueue, String toolName) {
        this.outputQueue = outputQueue;
        TOOL_NAME = toolName;
        FILE_NAME = fileName;
    }

    public FileWriterConsumer(String fileName, BlockingQueue<ArrayList<String>> outputQueue, String recordName, String toolName) {
        this.outputQueue = outputQueue;
        this.RECORD_NAME = recordName;
        TOOL_NAME = toolName;
        FILE_NAME = fileName;
    }

    @Override
    public void run() {
//        try {
////            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(java.io.FileDescriptor.out), "ASCII"), 512);
//            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(FILE_NAME))), BUFFER_SIZE);
//            ArrayList<String> list;
//            long outputCount = 0L;
//            while (!(list = outputQueue.take()).isEmpty()) {
//                outputCount += list.size();
//                for (String line : list) {
//                    out.write(line);
//                    out.newLine();
//                }
//            }
//            out.flush();
//            //output remianing stored records
//            Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(outputCount) + " " + RECORD_NAME + " written-out to "+FILE_NAME,TOOL_NAME);
//        } catch (InterruptedException | UnsupportedEncodingException e) {
//            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
//        } catch (IOException e) {
//            Reporter.report("[ERROR]", e.getMessage(), TOOL_NAME);
//        }

        String newline = System.lineSeparator(); //.getProperty("line.separator");
        BufferedWriter writer = null;
        BufferedWriter writer2 = null;
        try {
//            if (task == Task.WRITE_FASTQ_PE) {
//            if (!outputQueue.peek().isEmpty()) {
//                GZIPOutputStream outStream = new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R1.fastq.gz"));
//                writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R1.fastq.gz")), "UTF-8"));
//                GZIPOutputStream outStream2 = new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R2.fastq.gz"));
//            writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R2.fastq.gz")), "UTF-8"));
//            } else if (task == Task.WRITE_FASTQ_SE) {
//                GZIPOutputStream outStream = new GZIPOutputStream(new FileOutputStream("tmp_SE.fastq.gz"));
//                writer = new BufferedWriter(new OutputStreamWriter(outStream, "UTF-8"));
//            } else if (task == Task.WRITE_FASTA_PE) {
////                GZIPOutputStream outStream = new GZIPOutputStream(new FileOutputStream("tmp_R1.fastq.gz"));
////                writer = new BufferedWriter(new OutputStreamWriter(outStream, "UTF-8"));
////                GZIPOutputStream outStream2 = new GZIPOutputStream(new FileOutputStream("tmp_R2.fastq.gz"));
////                writer2 = new BufferedWriter(new OutputStreamWriter(outStream2, "UTF-8"));
//            } else if (task == Task.WRITE_FASTA_SE) {
////                GZIPOutputStream outStream = new GZIPOutputStream(new FileOutputStream("tmp_SE.fastq.gz"));
////                writer = new BufferedWriter(new OutputStreamWriter(outStream, "UTF-8"));
//            }

            ArrayList<String> list;
            long outputCount = 0L;
            while (!(list = outputQueue.take()).isEmpty()) {
                if (writer == null) {
                    writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R1.fastq.gz")), "UTF-8"));
                    writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(FILE_NAME + "_R2.fastq.gz")), "UTF-8"));
                }
                outputCount += list.size();
//                if (task == Task.WRITE_FASTQ_SE) {
//                    writer.write(read.replaceAll("\t", newline));
//                } else if (task == Task.WRITE_FASTQ_PE) {
                for (String read : list) {
                    StringBuilder sb = new StringBuilder();
                    StringBuilder sb2 = new StringBuilder();
                    String[] splits = read.split("\t");
                    if (splits.length == 8) {
                        for (int i = 0; i < 4; i++) {
                            sb.append(splits[i]).append(newline);
                            sb2.append(splits[i + 4]).append(newline);
                        }
                    } else {
                        System.err.println("Output writer error, expecting fastq record to have 4 or 8 tab delimited fields! Offending line:\n" + read);
                        System.err.println("Terminating!");
                        System.exit(1);
                    }
                    writer.write(sb.toString()); //                    writer.newLine();
                    writer2.write(sb2.toString());
                }
                writer.flush();
                writer2.flush();
            }
            if (outputCount > 0) {
                Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(outputCount) + " " + RECORD_NAME + " written-out to " + FILE_NAME, TOOL_NAME);
            }
//            }
//            queue.put("TERMINATE"); //inform other threads -should not be necessary - one thread for writing
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (writer != null) {
                    writer.close();
                }
                if (writer2 != null) {
                    writer2.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

    }

}
