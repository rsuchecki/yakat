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

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Orf implements Comparable<Orf> {

    private String parenId;
    private int from;
    private int to;
    private int frame;
    private boolean hasStopCodon;

    public Orf(String parenId, int from, int to, int frame, boolean hasStopCodon) {
        this.parenId = parenId;
        this.from = from;
        this.to = to;
        this.frame = frame;
        this.hasStopCodon = hasStopCodon;
    }

    public String getParenId() {
        return parenId;
    }

    public long getFrom() {
        return from;
    }

    public long getTo() {
        return to;
    }

    public int getFrame() {
        return frame;
    }

    public long getLength() {
        return getTo() - getFrom() + 1;
    }

    public boolean hasStopCodon() {
        return hasStopCodon;
    }


    @Override
    public int compareTo(Orf o) {
        return (int) (o.getLength() - getLength());
    }
    
    
    

}
