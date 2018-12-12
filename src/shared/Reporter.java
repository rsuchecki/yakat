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
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

/**
 * Generic reporting interface
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Reporter {

    private static int PRINT_WIDTH = 160;

    /**
     * Prints to std err, a message formatted as per example: 2015-06-24 Wed 10:10:57 [KmerExtender.jar] [INFO]
     * Initialized [1.12 MB used]
     *
     * @param level : use e.g. "[WARNING]", "[INFO]" etc.
     * @param message
     */
    public static void report(String level, String message, String tool) {
        System.err.print(formatReport(level, message, true, tool));
//        Date date = new Date(System.currentTimeMillis());
//        String newline = System.lineSeparator();
//        Runtime runtime = Runtime.getRuntime();
//        StringBuilder mem = new StringBuilder();
//        mem.append("JVM mem: ");
//        mem.append(CommonMaths.getBytesMultiple(runtime.totalMemory())).append(" total, ");
//        mem.append(CommonMaths.getBytesMultiple(runtime.freeMemory())).append(" free.");
//        long usedMem;
//        long total = runtime.totalMemory();
//        long free = runtime.freeMemory();
//        usedMem = total - free;
//        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd EEE HH:mm:ss");
//        String format = "%s %-20s %-10s %s %s" + newline;
//        
////        System.err.printf(format, dateFormat.format(date), "[KmerExtender.jar]", level, message, "[" + CommonMaths.getBytesMultiple(usedMem) + " used]");
////            System.err.printf(format, dateFormat.format(date), "KmerExtender.jar", level, message);
////            System.err.printf("[%s,\t%s]" + "\t%s\t%s" + newline, dateFormat.format(date), mem.toString(), message, CommonMaths.getBytesMultiple(usedMem));
    }

    /**
     * Prints to std err, a message formatted as per example: 2015-06-24 Wed 10:10:57 [KmerExtender.jar] [INFO]
     * Initialized
     *
     * @param level : use e.g. "[WARNING]", "[INFO]" etc.
     * @param message
     * @param tool
     */
    public static void reportNoMem(String level, String message, String tool) {
        System.err.print(formatReport(level, message, false, tool));
    }

    public static String formatReport(String level, String message, String tool) {
        return formatReport(level, message, true, tool);
    }

    public static String formatReport(String level, String message, boolean printMemoryUsage, String tool) {
        Date date = new Date(System.currentTimeMillis());
        String newline = System.lineSeparator();
        Runtime runtime = Runtime.getRuntime();
        StringBuilder mem = new StringBuilder();
        mem.append("JVM mem: ");
        mem.append(CommonMaths.getBytesMultiple(runtime.totalMemory())).append(" total, ");
        mem.append(CommonMaths.getBytesMultiple(runtime.freeMemory())).append(" free.");
        long usedMem;
        long total = runtime.totalMemory();
        long free = runtime.freeMemory();
        usedMem = total - free;
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd EEE HH:mm:ss");

        String terminalColumns = System.getenv("COLUMNS");
        if (terminalColumns != null) {
            PRINT_WIDTH = Integer.parseInt(terminalColumns);
        }

//        int helpLineWidth = PRINT_WIDTH - offset;
        String format = "%s %-22s %-10s %s %s" + newline;        
        String memUse = "";
        if (printMemoryUsage) {
            memUse = "[" + CommonMaths.getBytesMultiple(usedMem) + "]";
        }

//        int prefixLen = dateFormat.format(date).length() + 1;
//        int suffixLen = 1 + memUse.length();
//        int coreLen = prefixLen + suffixLen;
//        int expectedLen = coreLen + message.length();
//        if (expectedLen > PRINT_WIDTH) {
//            message = wrapString(newline+message, PRINT_WIDTH - prefixLen, prefixLen, 33);
////            message += newline;
////            message = wrapString(message, PRINT_WIDTH - coreLen, prefixLen, memUse);
//        }

        String s = String.format(format, dateFormat.format(date), "[" + tool + "]", level, message, memUse);

//        System.err.println("len="+s.length());
//        System.err.println("wdt="+PRINT_WIDTH);
//        System.err.println(s);
//        System.exit(1);
        return s;
    }

    public static void writeToFile(String fName, String contents, boolean append) {
        if (fName != null) {
            try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName, append)))) {
                out.println(contents);
            } catch (IOException e) {
                report("[ERROR]", e.getMessage(), "");
            }
        }
    }

    public static void writeToFile(String fName, ArrayList<String> contents, boolean append) {
        if (fName != null) {
            try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName, append)))) {
                for (String line : contents) {
                    out.println(line);
                }
            } catch (IOException e) {
                report("[ERROR]", e.getMessage(), "");
            }
        }
    }

    public static void writeHistogramToFile(String fName, int[] contents, boolean append, int k) {
        if (fName != null) {
            try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName, append)))) {
                out.println("length\tfrequency");
                for (int i = 1; i < contents.length; i++) {
                    if (i > k) {
                        out.println(i + "\t" + contents[i]);
                    }
                }
                out.println(">max\t" + contents[0]);

            } catch (IOException e) {
                report("[ERROR]", e.getMessage(), "");
            }
        }
    }

    public static String wrapString(String s, int lineLength) {
        StringBuilder sb = new StringBuilder();
        String[] toks = s.split(" |\t|\n");
        int length = 0;
        for (String tok : toks) {
            if (length + tok.length() > lineLength) {
                sb.append(System.lineSeparator());
                length = 0;
            }
            sb.append(tok).append(" ");
            length += tok.length() + 1;
        }
        return sb.toString();
    }

    public static String wrapString(String s, int lineLength, int offset) {
        String offsetString = String.format("%" + offset + "s", " ");
        StringBuilder sb = new StringBuilder();
        String[] toks = s.split(" |\t|\n");
        int length = 0;
        for (String tok : toks) {
            if (length + tok.length() > lineLength - offset) {
                sb.append(System.lineSeparator());
                sb.append(offsetString);
                length = 0;
            }
            sb.append(tok).append(" ");
            length += tok.length() + 1;
        }
        return sb.toString();
    }
    
    public static String wrapString(String s, int lineLength, int offset, int extraOffsetFirstLine) {
        String offsetString = String.format("%" + offset + "s", " ");
        StringBuilder sb = new StringBuilder();
        String[] toks = s.split(" |\t|\n");
        int length = 0;
        for (String tok : toks) {
            if (length + tok.length() > lineLength - offset - extraOffsetFirstLine) {
                extraOffsetFirstLine = 0;
                sb.append(System.lineSeparator());
                sb.append(offsetString);
                length = 0;
            }
            sb.append(tok).append(" ");
            length += tok.length() + 1;
        }
        return sb.toString();
    }
    
    public static String wrapString(String s, int lineLength, int offset, String appendToFirstLine) {
        String offsetString = String.format("%" + offset + "s", " ");
        StringBuilder sb = new StringBuilder();
        String[] toks = s.split(" |\t|\n");
        int length = 0;
        int addToFirst = appendToFirstLine.length()+1;
        for (String tok : toks) {
            if (length + tok.length() > lineLength - offset - addToFirst) {                
                sb.append(appendToFirstLine);
                appendToFirstLine="";
                addToFirst = 0;
                sb.append(System.lineSeparator());
                sb.append(offsetString);
                length = 0;
            }
            sb.append(tok).append(" ");
            length += tok.length() + 1;
        }
        return sb.toString();
    }
}
