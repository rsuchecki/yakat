/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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

import java.util.regex.Pattern;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Codon {    
    private final String seq;
    private final int pos;
    private final int frame;
    
    
    public Codon(String seq, int pos, int frame) {
        this.seq = seq;
        this.pos = pos;
        this.frame = frame;
    }

    public String getSeq() {
        return seq;
    }

    public int getPosition() {
        return pos;
    }

    public int getFrame() {
        return frame;
    }
    
    public boolean isStart() {
        return (frame > 0 && getSeq().equalsIgnoreCase("ATG")) || (frame < 0 && getSeq().equalsIgnoreCase("CAT"));
    }
    
    public boolean isStop(Pattern stopForward, Pattern stopReverse) {
//        Pattern p = Pattern.compile("");
//        p.matcher(getSeq().toUpperCase()).;
//        return (frame > 0 && getSeq().toUpperCase().matches("T(AA|AG|GA)")) || (frame < 0 && getSeq().toUpperCase().matches("(TT|TC|CT)A"));
        return (frame > 0 && stopForward.matcher(getSeq()).matches()) || (frame < 0 && stopReverse.matcher(getSeq()).matches());
    }
}
