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
import java.util.concurrent.atomic.AtomicInteger;
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
    private final String DIR_NAME;
    private final String FILE_NAME;
    private int PRODUCER_THREADS;
    private final String R1_SUFFIX;
    private final String R2_SUFFIX;
    private final String SE_SUFFIX;


    public FileWriterConsumer(BlockingQueue<ArrayList<String>> outputQueue, String TOOL_NAME, String DIR_NAME, String FILE_NAME, int PRODUCER_THREADS, String R1_SUFFIX, String R2_SUFFIX, String SE_SUFFIX) {
        this.outputQueue = outputQueue;
        this.TOOL_NAME = TOOL_NAME;
        this.DIR_NAME = DIR_NAME;
        this.FILE_NAME = FILE_NAME;
        this.PRODUCER_THREADS = PRODUCER_THREADS;
        this.R1_SUFFIX = R1_SUFFIX;
        this.R2_SUFFIX = R2_SUFFIX;
        this.SE_SUFFIX = SE_SUFFIX;
    }
    
    

//    public FileWriterConsumer(String fileName, BlockingQueue<ArrayList<String>> outputQueue, String recordName, String toolName) {
//        this.outputQueue = outputQueue;
//        this.RECORD_NAME = recordName;
//        TOOL_NAME = toolName;
//        FILE_NAME = fileName;
//    }
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
        BufferedWriter writer1 = null;
        BufferedWriter writer2 = null;
        BufferedWriter writerOrphans = null;
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
            
            while (!(list = outputQueue.take()).isEmpty() || --PRODUCER_THREADS > 0) {

                if (writer1 == null && !list.isEmpty()) {
                    //create 2 writers                    
                    writer1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(DIR_NAME+ File.separator +FILE_NAME+R1_SUFFIX)), "UTF-8"), BUFFER_SIZE);
                    writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(DIR_NAME+ File.separator +FILE_NAME+R2_SUFFIX)), "UTF-8"), BUFFER_SIZE);
                }

//                    PRODUCER_THREADS.decrementAndGet();
//                    continue;
//                }
                outputCount += list.size();
//                if (task == Task.WRITE_FASTQ_SE) {
//                    writer.write(read.replaceAll("\t", newline));
//                } else if (task == Task.WRITE_FASTQ_PE) {
                while (!list.isEmpty()) {
                    String line = list.remove(list.size() - 1);
                    StringBuilder sb1 = new StringBuilder();
                    StringBuilder sb2 = new StringBuilder();
                    StringBuilder sbOrphans = new StringBuilder();
                    String[] splits = line.split("\t");
                    if (splits.length == 4) {
                        for (int i = 0; i < 4; i++) {
                            sbOrphans.append(splits[i]).append(newline);
                        }
                        if (writerOrphans == null) {
                            writerOrphans = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(DIR_NAME+ File.separator +FILE_NAME+SE_SUFFIX)), "UTF-8"), BUFFER_SIZE);
                        }
                        writerOrphans.write(sbOrphans.toString());
                    } else if (splits.length == 8) {
                        for (int i = 0; i < 4; i++) {
                            sb1.append(splits[i]).append(newline);
                            sb2.append(splits[i + 4]).append(newline);
                        }
                        writer1.write(sb1.toString());
                        writer2.write(sb2.toString());
                    } else {
                        System.err.println("Output writer error, expecting fastq record to have 4 or 8 tab delimited fields! Offending line:\n" + line);
                        System.err.println("Terminating!");
                        System.exit(1);
                    }
                }
                if (writer1 != null) {
                    writer1.flush();
                    writer2.flush();
                }
                if (writerOrphans != null) {
                    writerOrphans.flush();
                }
            }
            if (outputCount > 0) {
                Reporter.report("[INFO]", NumberFormat.getNumberInstance().format(outputCount) + " " + RECORD_NAME + " written-out to " + FILE_NAME, TOOL_NAME);
//                Reporter.report("[INFO]", "[Thread " + Thread.currentThread().getId()+"] "+NumberFormat.getNumberInstance().format(outputCount) + " " + RECORD_NAME + " written-out to " + FILE_NAME, TOOL_NAME);
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (writer1 != null) {
                    writer1.close();
                }
                if (writer2 != null) {
                    writer2.close();
                }
                if (writerOrphans != null) {
                    writerOrphans.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

    }

}
