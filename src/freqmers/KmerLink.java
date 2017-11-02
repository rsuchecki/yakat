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
package freqmers;

import shared.Sequence;

/**
 *
 * @author rad
 */
public class KmerLink implements Comparable<KmerLink>{    
    private final KmerFilter kmerFilter;
    private final int startPosition;
    private boolean unique = true;

    public KmerLink(KmerFilter snpFilter, int startPosition) {
        this.kmerFilter = snpFilter;
        this.startPosition = startPosition;
    }

    public KmerFilter getSnpFilter() {
        return kmerFilter;
    }
    
    public int getStartPosition() {
        return startPosition;
    }
    
    public Sequence getParentSequence() {
        return kmerFilter.getSequence();
    }
    
    public boolean setMer(short value) {
            return kmerFilter.setMer(getStartPosition(), value);
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
