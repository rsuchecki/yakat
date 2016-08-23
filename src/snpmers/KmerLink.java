/*
 * Copyright 2016 rad.
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
package snpmers;

import shared.Sequence;

/**
 *
 * @author rad
 */
public class KmerLink implements Comparable<KmerLink>{    
    private final SnpFilter snpFilter;
    private final boolean parentOne;        
    private final int startPosition;
    private final boolean reverseComplement;
    private boolean unique = true;

    public KmerLink(SnpFilter snpFilter, boolean parentOne, int startPosition, boolean reverseComplement) {
        this.snpFilter = snpFilter;
        this.parentOne = parentOne;
        this.startPosition = startPosition;
        this.reverseComplement = reverseComplement;
    }

    public SnpFilter getSnpFilter() {
        return snpFilter;
    }
    
    public int getStartPosition() {
        return startPosition;
    }
    
    public Sequence getParentSequence() {
        return parentOne ? snpFilter.getSequence1() : snpFilter.getSequence2();
    }
    
    public boolean setMer(short value) {
        if(parentOne) {
            return snpFilter.setMer1(getStartPosition(), value);
        } else {
            return snpFilter.setMer2(getStartPosition(), value);            
        }        
    }

    public boolean isParentOne() {
        return parentOne;
    }

    
    public boolean isReverseComplement() {
        return reverseComplement;
    }

    
    
    public boolean isUnique() {
        return unique;
    }

    public void setUnique(boolean unique) {
        this.unique = unique;
    }
    

    @Override
    public int compareTo(KmerLink o) {
        return getStartPosition() - o.getStartPosition();
    }
    
    
    
}
