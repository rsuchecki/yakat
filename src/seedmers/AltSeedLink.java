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
package seedmers;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class AltSeedLink {
    private final Seed parent;
    private final int position;
    private final byte base;
    private final boolean forward;
    
    public AltSeedLink(Seed parent, int position, char base, boolean forward) {
        this.parent = parent;
        this.position = position;
        this.base = encodeChar(base);
        this.forward = forward;
    }

    public int getPosition() {
        return position;
    }

    public Seed getParent() {
        return parent;
    }
    
    private byte encodeChar(char base) {
        switch(base) {
            case 'A' : return 0;
            case 'C' : return 1;
            case 'G' : return 2;
            case 'T' : return 3;
        }
        return -1;
    }

    public byte getBase() {
        return base;
    }
    
    
    public  char getBaseChar() {
        switch(base) {
            case 0 : return 'A';
            case 1 : return 'C';
            case 2 : return 'G';
            case 3 : return 'T';
        }
        return '#';
    }

    public boolean isForward() {
        return forward;
    }
    
    
    
}
