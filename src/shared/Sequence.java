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

    public ArrayList<Orf> getOrfs() {
        ArrayList<Codon> codons = new ArrayList<>();
        Pattern stopCodons = Pattern.compile("ATG|(T(AA|AG|GA))", Pattern.CASE_INSENSITIVE);
        int frame = 0;
        for (int i = 0; i < getLength() - 2; i++) {
            frame = frame < 3 ? frame + 1 : 1;
            String codon = getSequenceString().substring(i, i + 3);
            Matcher matcher = stopCodons.matcher(codon);
            if (matcher.matches()) {
                codons.add(new Codon(codon, (i + 1), frame));
            }
        }
        //accomodating for incomplete ORFs - missing a stop codon        
        for (frame = 1; frame <= 3; frame++) {
            codons.add(new Codon("last", getLength(), frame));
        }
        //NOW RC 
        stopCodons = Pattern.compile("CAT|((TT|TC|CT)A)", Pattern.CASE_INSENSITIVE);
        frame = 0;
        for (int i = getLength() - 3; i >= 0; --i) {
            frame = frame > -3 ? frame - 1 : -1;
            String codon = getSequenceString().substring(i, i + 3);
            Matcher matcher = stopCodons.matcher(codon);
            if (matcher.matches()) {
                codons.add(new Codon(codon, (i + 1), frame));
            }
        }
        //accomodating for incomplete ORFs - missing a stop codon        
        for (frame = -3; frame <= -1; frame++) {
            codons.add(new Codon("last", 1, frame));
        }
        ArrayList<Orf> orfs = new ArrayList<>();
        HashMap<Integer, Codon> previousCodonInFrame = new HashMap<>();
        for (Codon codon : codons) {
//            System.err.println(getId() + "\t" + codon.isStart() + "\t" + codon.getPosition() + "\t" + codon.getFrame() + "\t" + codon.getSeq());
            int currentFrame = codon.getFrame();
            if (previousCodonInFrame.containsKey(codon.getFrame())) {
                Codon previous = previousCodonInFrame.get(codon.getFrame());
                if (previous.isStart()) {
                    if (codon.isStop() || codon.getSeq().equalsIgnoreCase("last")) {
                        int left = currentFrame > 0 ? previous.getPosition() : codon.getPosition();
                        int right = currentFrame > 0 ? codon.getPosition() + 2 : previous.getPosition() + 2;
//      print $1,($2=="last" && $4<0)?"<"left:left,($2=="last" && $4>0)?">"right:right,$4,right-left+1 
                        orfs.add(new Orf(getId(), left, right, currentFrame, codon.isStop()));
                        previousCodonInFrame.put(currentFrame, codon);
                    } else if (codon.isStart()) {
                        //ignore start codon contained in ORF
                        continue;
                    }
                }
            }
            previousCodonInFrame.put(currentFrame, codon);
        }
        return orfs;
    }
}
