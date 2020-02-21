/*
 * Copyright 2017 Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au.
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
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class FastaIndexed {

    //FAI index fields:
    //NAME          Name of this reference sequence
    //LENGTH	Total length of this reference sequence, in bases
    //OFFSET	Offset within the FASTA file of this sequence's first base
    //LINEBASES	The number of bases on each line
    //LINEWIDTH	The number of bytes in each line, including the newline            
    private String TOOL_NAME;
    private String fastaFile;
    private HashMap<String, Integer> lengthsMap = new HashMap<>();
    private HashMap<String, Integer> offsetsMap = new HashMap<>();
    private HashMap<String, Integer> lineBasesMap = new HashMap<>();
    private HashMap<String, Integer> lineWidthsMap = new HashMap<>();
    private ArrayList<String> ids;
    
    public FastaIndexed(String TOOL_NAME, String fasta, String fai) {
        this.ids = new ArrayList<>();
        this.TOOL_NAME = TOOL_NAME;
        fastaFile = fasta;
        if (fai == null) {
            //TODO create index
        } else {
            readIndex(fai);
        }
    }

    
    private void readIndex(String fai) {
        BufferedReader index = null;
        try {
            index = new BufferedReader(new FileReader(new File(fai)));
            String line;
            while ((line = index.readLine()) != null && !line.isEmpty()) {
                String[] toks = line.split("\t");
                ids.add(toks[0]);
                lengthsMap.put(toks[0], Integer.parseInt(toks[1]));
                offsetsMap.put(toks[0], Integer.parseInt(toks[2]));
                lineBasesMap.put(toks[0], Integer.parseInt(toks[3]));
                lineWidthsMap.put(toks[0], Integer.parseInt(toks[4]));
//                System.out.println(line);
//                if(!Objects.equals(lineBasesMap.get(toks[0]), lengthsMap.get(toks[0]))) {
//                    Reporter.report("[ERROR]", "Wrapped FASTA not supported. Offending record: "+line , TOOL_NAME);
//                    System.exit(1);
//                }              
            }
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found exception: " + ex.getMessage(), TOOL_NAME);
            System.exit(1);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        } finally {
            try {
                if (index != null) {
                    index.close();
                }
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }

        }

    }
    
    

    public ArrayList<String> getIds() {
        return ids;
    }


    
    
    public String getSequence(String id) {
        return getSequence(id, null, null);
    }

    /**
     * Extract
     *
     * @param id
     * @param from
     * @param to
     * @return
     */
    public String getSequence(String id, Integer from, Integer to) {
        RandomAccessFile file = null;
        byte[] buff = null;
        try {
            file = new RandomAccessFile(fastaFile, "r");
            Integer start = offsetsMap.get(id);
            if(start == null) {
                Reporter.report("[ERROR]", id+" <- identifier  not found in .fai index for: "+fastaFile, TOOL_NAME);                
                System.exit(1);
            }
            if(from != null) {
                start += from - 1;
            }
            file.seek(start);
           
           
            int len = to == null ? lengthsMap.get(id) : to - from + 1;               
            buff = new byte[len];
            
            //Wrapped FASTA
            int lineWidth = lineWidthsMap.get(id);
            int lineBases = lineBasesMap.get(id);
            int wrapBytes = lineWidth-lineBases;
//            byte[] sink = new byte[wrapBytes];
            if(lineWidth < len) {
                int lines = len / lineWidth;
//                len += lines;
                int offset = 0;
                for (int i = 0; i < lines; i++) {                     
                   offset += file.read(buff,offset,lineBases);
//                   byte tmp[] = buff.clone();
//                   String sofar = new String(tmp);
//                   file.read(sink);
                   file.skipBytes(wrapBytes);
                }
                file.read(buff,offset,buff.length-offset);

//                    offset = file.read(buff,0,lineBases);
//
//                      file.skipBytes(wrapBytes);
////                    byte newlineOrReturn = file.readByte();
////                    if(newlineOrReturn != 10) { //If not \n
////                        if(newlineOrReturn == 13) { //If \r
////                            newlineOrReturn = file.readByte();
////                        }
////                    }
////                    String should_be_newline = Byte.toString(file.readByte());
//                    offset += file.read(buff,offset,lineBases);                    
//                    file.skipBytes(wrapBytes);
//                    offset+= file.read(buff,offset,lineBases);
//                    
//                    file.skipBytes(wrapBytes);
//                    file.read(buff,offset,buff.length-offset);
                    
//                }
//                String read = new String(buff);
//                int  readLen = read.length();
//                
//                System.err.println(read);
//                System.err.println(readLen);
                
//                for (int i = 0; i < buff.length; i++) {
//                    System.out.printf("%4d", i);
//                }
//                for (int i = 0; i < buff.length; i++) {
//                    System.out.printf("%4s", buff[i]);
//                }
//                System.out.println();
//                int x =0;
            } else {                       
                
                file.read(buff);
            }
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found exception: " + ex.getMessage(), TOOL_NAME);
            System.exit(1);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        } finally {
            try {
                if (file != null) {
                    file.close();
                }
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
        return buff == null ? "" : new String(buff);
//        return buff == null ? "" : (Objects.equals(lineBasesMap.get(id), lengthsMap.get(id)) ? new String(buff) : new String(buff).replaceAll("\\r\\n?|\\n", ""));
    }
}
