/*
 * Copyright 2018 rad.
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
 * Tobe stored in auxiliary map for a PairMer which has alternative clip(s)
 * @author rad suchecki
 */
public class PairMerCLips {

    private byte clipCountLeftA;
    private byte clipCountLeftC;
    private byte clipCountLeftG;
    private byte clipCountLeftT;
    private byte clipCountRightA;
    private byte clipCountRightC;
    private byte clipCountRightG;
    private byte clipCountRightT;

    public void increment(byte base, boolean left, byte freq) {

        switch(base) {
            case 8 : {                
                if(left) {
                    setClipCountLeftA(increment(getClipCountLeftA(), freq));
                } else {
                    setClipCountRightA(increment(getClipCountRightA(), freq));                    
                }
            }
            case 4 : {                
                if(left) {
                    setClipCountLeftC(increment(getClipCountLeftC(), freq));
                } else {
                    setClipCountRightC(increment(getClipCountRightC(), freq));                    
                }
            }
            case 2 : {                
                if(left) {
                    setClipCountLeftG(increment(getClipCountLeftG(), freq));
                } else {
                    setClipCountRightG(increment(getClipCountRightG(), freq));                    
                }
            }
            case 1 : {                
                if(left) {
                    setClipCountLeftT(increment(getClipCountLeftT(), freq));
                } else {
                    setClipCountRightT(increment(getClipCountRightT(), freq));                    
                }
            }
            default: Reporter.reportNoMem("[ERROR]", "Unexpected, non-ACGT clip value ", this.getClass().getCanonicalName());
        }
        
        
        
    }

    public byte getClipCountLeftA() {
        return clipCountLeftA;
    }

    public byte getClipCountLeftC() {
        return clipCountLeftC;
    }

    public byte getClipCountLeftG() {
        return clipCountLeftG;
    }

    public byte getClipCountLeftT() {
        return clipCountLeftT;
    }

    public byte getClipCountRightA() {
        return clipCountRightA;
    }

    public byte getClipCountRightC() {
        return clipCountRightC;
    }

    public byte getClipCountRightG() {
        return clipCountRightG;
    }

    public byte getClipCountRightT() {
        return clipCountRightT;
    }

    public void setClipCountLeftA(byte clipCountLeftA) {
        this.clipCountLeftA = clipCountLeftA;
    }

    public void setClipCountLeftC(byte clipCountLeftC) {
        this.clipCountLeftC = clipCountLeftC;
    }

    public void setClipCountLeftG(byte clipCountLeftG) {
        this.clipCountLeftG = clipCountLeftG;
    }

    public void setClipCountLeftT(byte clipCountLeftT) {
        this.clipCountLeftT = clipCountLeftT;
    }

    public void setClipCountRightA(byte clipCountRightA) {
        this.clipCountRightA = clipCountRightA;
    }

    public void setClipCountRightC(byte clipCountRightC) {
        this.clipCountRightC = clipCountRightC;
    }

    public void setClipCountRightG(byte clipCountRightG) {
        this.clipCountRightG = clipCountRightG;
    }

    public void setClipCountRightT(byte clipCountRightT) {
        this.clipCountRightT = clipCountRightT;
    }
    
    
    
    private byte increment(byte previous, byte value) {
        return (byte) Math.min(previous + value, Byte.MAX_VALUE);
    }
}
