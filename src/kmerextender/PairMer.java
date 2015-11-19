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
package kmerextender;

import shared.Reporter;

/**
 * Superclass - not to be used directly holds up to 2 k-mers paired
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMer {//implements Comparable<PairMer> {

    //Object overhead 8 B
//    private int[] kmerCoreBitsArray;  //12B + len*4 Bytes 
    private char clipLeft = '#';   //2B                     ///*TODO*/: encode to 2-3b if sticking to int array
    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
    private byte storedCount;  //1B
    private boolean invalid;  //1B
    private boolean visited;//1B
    //then round to multi of 8

    //could save 7B by encoding all the fields in a single B (storedCount would not be stored exactly but the state of the other bits could indicate is count is 1,2 or more)
    //7 (sign) 
    //6,5,4    left clip:  \fi,A,C,G,T or X for extensions conflict
    //3,2,1    right clip: \fi,A,C,G,T or X for extensions conflict
    //0        wasVisited: 0==false, 1==true 
    
    
    /**
     * It is assumed that the input PairMer matches this one (another.core ==
     * this.core)
     *
     * @param another : newly generated PairMer holding a single k-mer
     * @param inputKmersUnique : set true if no duplicate k-mers expected
     */
    public synchronized void addKmerSynchronized(PairMer another, boolean inputKmersUnique) {
        if (isInvalid() || (inputKmersUnique && getStoredCount() > 1) || (!inputKmersUnique && getStoredCount() > 2)) { //if already invalid PairMer  or more than second kmer being added
            setIsInvalid();
        } else {
            //If input k-mers are non-unique, ie, a k-mer may appear more than once, we need to ensure that we ignore it            
            if (!inputKmersUnique) {
                if ((hasLeftClip() && getClipLeft() == another.getClipLeft()) || (hasRightClip() && getClipRight() == another.getClipRight())) {
                    //fine, new k-mer is identical to a k-mer already stored
//                    System.err.println("Adding same kmer again");
                } else { //it is a different k-mer
                    if ((hasLeftClip() && another.hasLeftClip()) || (hasRightClip() && another.hasRightClip())) {
                        setIsInvalid();
                    } else if (!hasLeftClip() && another.hasLeftClip()) {
                        setClipLeft(another.getClipLeft());
                    } else if (!hasRightClip() && another.hasRightClip()) {
                        setClipRight(another.getClipRight());
                    } else {
                        Reporter.report("[BUG?]", "Unexpected [2], addKmer() at " , this.getClass().getSimpleName());
                    }
                    incrementStoredCount();
                }
            } else { //i.e. input k-mers are unique (no duplicates which we need to ignore)
                if (hasBothClips()) {
                    setIsInvalid();
                } else if ((hasLeftClip() && another.hasLeftClip()) || (hasRightClip() && another.hasRightClip())) {
                    setIsInvalid();
                } else if (!hasLeftClip() && another.hasLeftClip()) {
                    setClipLeft(another.getClipLeft());
                } else if (!hasRightClip() && another.hasRightClip()) {
                    setClipRight(another.getClipRight());
                } else {
                    Reporter.report("[BUG?]", "Unexpected [3], addKmer() at " , this.getClass().getSimpleName());
                }
                incrementStoredCount();
            }
        }
    }

    protected boolean hasLeftClip() {
        return clipLeft != '#';
    }

    protected boolean hasRightClip() {
        return clipRight != '#';
    }

    protected boolean hasBothClips() {
        return hasLeftClip() && hasRightClip();
    }

    protected void incrementStoredCount() {
        if (storedCount < Byte.MAX_VALUE) {
            storedCount++;
        }
    }

    protected char getClipLeft() {
        return clipLeft;
    }

    protected char getClipRight() {
        return clipRight;
    }

    /**
     * Placeholder method implemented by subclasses
     *
     * @param coreLength
     * @return
     */
    public String decodeCore(int coreLength) {
        return null;
    }

    public String getPairMerString(int k) {
        StringBuilder sb = new StringBuilder();
        sb.append(getClipLeft());
        sb.append(decodeCore(k - 1));
        sb.append(getClipRight());
        return sb.toString();
    }

    public boolean isInvalid() {
        return invalid;
    }

    protected byte getStoredCount() {
        return storedCount;
    }

    protected void setClipLeft(char clipLeft) {
        this.clipLeft = clipLeft;
    }

    protected void setClipRight(char clipRight) {
        this.clipRight = clipRight;
    }

    protected void setIsInvalid() {
        this.invalid = true;
    }

    public boolean isVisited() {
        return visited;
    }

    public void setVisited() {
        this.visited = true;
    }

}