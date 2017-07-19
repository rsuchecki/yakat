/*
 * Copyright 2016 Australian Centre For Plant Functional Genomics
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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Sequence {

    private String id;
    private String sequenceString;

    public Sequence(String id, String sequenceString) {
        this.id = id;
        this.sequenceString = sequenceString;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSequenceString() {
        return sequenceString;
    }

    public String getUnpaddedSequenceString() {
        return sequenceString.replaceAll("-", "");
    }

    public void setSequenceString(String sequenceString) {
        this.sequenceString = sequenceString;
    }

    public int getLength() {
        return sequenceString.length();
    }

    public int getLengthUnpadded() {
        return getUnpaddedSequenceString().length();
    }

    public CharSequence getFasta(boolean addLengthInDescriptionLine) {
        StringBuilder sb = new StringBuilder(">");
        sb.append(getId());
        if (addLengthInDescriptionLine) {
            sb.append(" ").append(getLength());
        }
        sb.append(System.lineSeparator()).append(getSequenceString());
        return sb;
    }

//    private void processCodon(Codon[] prevCodonIn, Codon currentCodon) {
//        int frame = currentCodon.getFrame();
//        int left = frame > 0 ? prevCodonIn[frame].getPosition() : currentCodon.getPosition();
//        int right = currentCodon.getFrame() > 0 ? currentCodon.getPosition() + 2 : prevCodonIn[frame].getPosition() + 2;
//        if (right - left + 1 >= minLen) {
//            orfs1.add(new Orf(getId(), left, right, frame, true));
//        }
//    }

    public ArrayList<Orf> getOrfs(int minLen, Pattern startStopCodonsForward, Pattern startStopCodonsReverse, Pattern stopForward, Pattern stopReverse) { //, Sequence parent, PrintStream bufferedOut) {
//        ArrayList<Codon> codons = new ArrayList<>(getLength() / 4); //ballpark size to prvent array copying on full
//        Pattern startStopCodons = Pattern.compile("ATG|(T(AA|AG|GA))", Pattern.CASE_INSENSITIVE);        
        ArrayList<Orf> orfs1 = new ArrayList<>(getLength() / 40); //ballpark size to prvent array copying on full

        int frame = 0;
        Codon[] prevCodonIn = new Codon[4];
        for (int i = 0; i < getLength() - 2; i++) {
            frame = frame < 3 ? frame + 1 : 1;
//            String codon = getSequenceString().substring(i, i + 3);
//            Matcher matcher = startStopCodonsForward.matcher(codon);
            Matcher matcher = startStopCodonsForward.matcher(getSequenceString().subSequence(i, i + 3));
            if (matcher.matches()) {
                //for 2-step approach
                Codon codon = new Codon(getSequenceString().substring(i, i + 3), (i + 1), frame);
//                codons.add(codon);

                if (prevCodonIn[frame] != null && prevCodonIn[frame].isStart()) {
                    if (codon.isStop(stopForward, stopReverse)) {
                        int left = frame > 0 ? prevCodonIn[frame].getPosition() : codon.getPosition();
                        int right = frame > 0 ? codon.getPosition() + 2 : prevCodonIn[frame].getPosition() + 2;
                        if (right - left + 1 >= minLen) {
                            orfs1.add(new Orf(getId(), left, right, frame, true));
                        }
//                        previousCodonInFrame.put(currentFrame, codon);
                    } else if (codon.isStart()) {
                        //ignore start codon contained in ORF
                        continue;
                    }
                }
                prevCodonIn[frame] = codon;
            }
        }
        //accomodating for incomplete ORFs - missing a stop codon        
        for (frame = 1; frame <= 3; frame++) {
//            codons.add(new Codon("last", getLength(), frame));
            if (prevCodonIn[frame] != null && prevCodonIn[frame].isStart()) {
                int left = frame > 0 ? prevCodonIn[frame].getPosition() : getLength();
                int right = frame > 0 ? getLength() + 2 : prevCodonIn[frame].getPosition() + 2;
                if (right - left + 1 >= minLen) {
                    orfs1.add(new Orf(getId(), left, right, frame, false));
                }
            }
        }
        //NOW RC, reset arr
        prevCodonIn = new Codon[4];
//        startStopCodons = Pattern.compile("CAT|((TT|TC|CT)A)", Pattern.CASE_INSENSITIVE);
        frame = 0;
        for (int i = getLength() - 3;
                i >= 0; --i) {
            frame = frame > -3 ? frame - 1 : -1;
//            String codon = getSequenceString().substring(i, i + 3);
            Matcher matcher = startStopCodonsReverse.matcher(getSequenceString().subSequence(i, i + 3));
            if (matcher.matches()) {
//                codons.add(new Codon(getSequenceString().substring(i, i + 3), (i + 1), frame));
                Codon codon = new Codon(getSequenceString().substring(i, i + 3), (i + 1), frame);
                if (prevCodonIn[-frame] != null && prevCodonIn[-frame].isStart()) {
                    if (codon.isStop(stopForward, stopReverse)) {
                        int left = frame > 0 ? prevCodonIn[-frame].getPosition() : codon.getPosition();
                        int right = frame > 0 ? codon.getPosition() + 2 : prevCodonIn[-frame].getPosition() + 2;
                        if (right - left + 1 >= minLen) {
                            orfs1.add(new Orf(getId(), left, right, frame, true));
                        }
//                        previousCodonInFrame.put(currentFrame, codon);
                    } else if (codon.isStart()) {
                        //ignore start codon contained in ORF
                        continue;
                    }
                }
                prevCodonIn[-frame] = codon;
            }
        }
        //accomodating for incomplete ORFs - missing a stop codon        
        for (frame = -3;
                frame <= -1; frame++) {
//            codons.add(new Codon("last", 1, frame));
            if (prevCodonIn[-frame] != null && prevCodonIn[-frame].isStart()) {
                int left = frame > 0 ? prevCodonIn[-frame].getPosition() : getLength();
                int right = frame > 0 ? getLength() + 2 : prevCodonIn[-frame].getPosition() + 2;
                if (right - left + 1 >= minLen) {
                    orfs1.add(new Orf(getId(), left, right, frame, false));
                }
            }
        }
        return orfs1;
    }

//    private void outputOrf(Sequence sequence, PrintStream bufferedOut, boolean translate) {
//        CharSequence seq = sequence.getSequenceString().subSequence(orf.getFrom() - 1, orf.getTo());
//        StringBuilder sb = new StringBuilder(orf.getFastaHeader());
//        sb.append(System.lineSeparator());
//        if (orf.getFrame() < 0) {
//            seq = SequenceOps.getReverseComplement(seq);
//        }
//        if (translate) {
//            sb.append(SequenceOps.translate(seq));
//        } else {
//            sb.append(seq);
//        }
//        bufferedOut.println(sb);
//    }
//    public ArrayList<Orf> getOrfsUpdated(int minLen, Pattern startStopCodonsForward, Pattern startStopCodonsReverse, Pattern stopForward, Pattern stopReverse) {
//        ArrayList<Codon> codons = new ArrayList<>(getLength() / 4); //ballpark size to prvent array copying on full
////        Pattern startStopCodons = Pattern.compile("ATG|(T(AA|AG|GA))", Pattern.CASE_INSENSITIVE);
//        int frame = 0;
//        for (int i = 0; i < getLength() - 2; i++) {
//            frame = frame < 3 ? frame + 1 : 1;
////            String codon = getSequenceString().substring(i, i + 3);
////            Matcher matcher = startStopCodonsForward.matcher(codon);
//            Matcher matcher = startStopCodonsForward.matcher(getSequenceString().subSequence(i, i + 3));
//            if (matcher.matches()) {
//                codons.add(new Codon(getSequenceString().substring(i, i + 3), (i + 1), frame));
//            }
//        }
//
//        //accomodating for incomplete ORFs - missing a stop codon        
//        for (frame = 1; frame <= 3; frame++) {
//            codons.add(new Codon("last", getLength(), frame));
//        }
//        //NOW RC 
////        startStopCodons = Pattern.compile("CAT|((TT|TC|CT)A)", Pattern.CASE_INSENSITIVE);
//        frame = 0;
//        for (int i = getLength() - 3; i >= 0; --i) {
//            frame = frame > -3 ? frame - 1 : -1;
////            String codon = getSequenceString().substring(i, i + 3);
//            Matcher matcher = startStopCodonsReverse.matcher(getSequenceString().subSequence(i, i + 3));
//            if (matcher.matches()) {
//                codons.add(new Codon(getSequenceString().substring(i, i + 3), (i + 1), frame));
//            }
//        }
//        //accomodating for incomplete ORFs - missing a stop codon        
//        for (frame = -3; frame <= -1; frame++) {
//            codons.add(new Codon("last", 1, frame));
//        }
////        System.err.println("total codons "+codons.size()+"/"+getLength()/4);
//        ArrayList<Orf> orfs = new ArrayList<>(codons.size() / 10); //ballpark size to prvent array copying on full
//        HashMap<Integer, Codon> previousCodonInFrame = new HashMap<>();
//        for (Codon codon : codons) {
////            System.err.println(getId() + "\t" + codon.isStart() + "\t" + codon.getPosition() + "\t" + codon.getFrame() + "\t" + codon.getSeq());
//            int currentFrame = codon.getFrame();
//            if (previousCodonInFrame.containsKey(codon.getFrame())) {
//                Codon previous = previousCodonInFrame.get(codon.getFrame());
//                if (previous.isStart()) {
//                    if (codon.isStop(stopForward, stopReverse) || codon.getSeq().equalsIgnoreCase("last")) {
//                        int left = currentFrame > 0 ? previous.getPosition() : codon.getPosition();
//                        int right = currentFrame > 0 ? codon.getPosition() + 2 : previous.getPosition() + 2;
////      print $1,($2=="last" && $4<0)?"<"left:left,($2=="last" && $4>0)?">"right:right,$4,right-left+1 
//                        if (right - left + 1 >= minLen) {
//                            orfs.add(new Orf(getId(), left, right, currentFrame, codon.isStop(stopForward, stopReverse)));
//                        }
//                        previousCodonInFrame.put(currentFrame, codon);
//                    } else if (codon.isStart()) {
//                        //ignore start codon contained in ORF
//                        continue;
//                    }
//                }
//            }
//            previousCodonInFrame.put(currentFrame, codon);
//        }
//        return orfs;
//    }
//    public ArrayList<Orf> getOrfsOLD(int minLen, Pattern startStopCodonsForward, Pattern startStopCodonsReverse, Pattern stopForward, Pattern stopReverse) {
//        ArrayList<Codon> codons = new ArrayList<>(getLength()/4); //ballpark size to prvent array copying on full
//        HashMap<Integer, Codon> previousCodonInFrame = new HashMap<>();
//        int frame = 0;
//        for (int i = 0; i < getLength() - 2; i++) {
//            frame = frame < 3 ? frame + 1 : 1;
//            String codon = getSequenceString().substring(i, i + 3);
//            Matcher matcher = startStopCodonsForward.matcher(codon);
//            if (matcher.matches()) {
//                codons.add(new Codon(codon, (i + 1), frame));
//                                
//            }
//        }
//        //accomodating for incomplete ORFs - missing a stop codon        
//        for (frame = 1; frame <= 3; frame++) {
//            codons.add(new Codon("last", getLength(), frame));
//        }
//        //NOW RC 
//        frame = 0;
//        for (int i = getLength() - 3; i >= 0; --i) {
//            frame = frame > -3 ? frame - 1 : -1;
//            String codon = getSequenceString().substring(i, i + 3);
//            Matcher matcher = startStopCodonsReverse.matcher(codon);
//            if (matcher.matches()) {
//                codons.add(new Codon(codon, (i + 1), frame));
//            }
//        }
//        //accomodating for incomplete ORFs - missing a stop codon        
//        for (frame = -3; frame <= -1; frame++) {
//            codons.add(new Codon("last", 1, frame));
//        }
////        System.err.println("total codons "+codons.size()+"/"+getLength()/4);
//        ArrayList<Orf> orfs = new ArrayList<>(codons.size()/10); //ballpark size to prvent array copying on full
//        for (Codon codon : codons) {
////            System.err.println(getId() + "\t" + codon.isStart() + "\t" + codon.getPosition() + "\t" + codon.getFrame() + "\t" + codon.getSeq());
//            int currentFrame = codon.getFrame();
//            if (previousCodonInFrame.containsKey(codon.getFrame())) {
//                Codon previous = previousCodonInFrame.get(codon.getFrame());
//                if (previous.isStart()) {
//                    if (codon.isStop(stopForward, stopReverse) || codon.getSeq().equalsIgnoreCase("last")) {
//                        int left = currentFrame > 0 ? previous.getPosition() : codon.getPosition();
//                        int right = currentFrame > 0 ? codon.getPosition() + 2 : previous.getPosition() + 2;
////      print $1,($2=="last" && $4<0)?"<"left:left,($2=="last" && $4>0)?">"right:right,$4,right-left+1 
//                        if (right - left + 1 >= minLen) {
//                            orfs.add(new Orf(getId(), left, right, currentFrame, codon.isStop(stopForward, stopReverse)));
//                        }
//                        previousCodonInFrame.put(currentFrame, codon);
//                    } else if (codon.isStart()) {
//                        //ignore start codon contained in ORF
//                        continue;
//                    }
//                }
//            }
//            previousCodonInFrame.put(currentFrame, codon);
//        }
//        return orfs;
//    }
}
