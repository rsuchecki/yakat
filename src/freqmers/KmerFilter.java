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

import java.util.HashMap;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class KmerFilter implements Comparable<KmerFilter>{
    private final String id;
    private final Sequence sequence;
    private short[] mers; //RECORD START POSITIONS OF ENCOUNTERED k-mers 
    private final HashMap<String, KmerFilterStats> samplesToStatsMap;
    private int mersCount;
    private int maxUniqMers;
    
    
    /**
     *
     * @param id
     * @param sequence1
     * @param k
     * @param TOOL_NAME
     */
    public KmerFilter(String id, Sequence sequence1, int k,  String TOOL_NAME) {
        this.sequence = sequence1;
        mers = new short[this.sequence.getLength()-k+1]; //
        this.id = id;
        samplesToStatsMap = new HashMap<>();
    }

    public void clearMers(String sample) {
        mers = new short[getSequence().getLength() + 1]; //
        mersCount = 0;
    }

    public Sequence getSequence() {
        return sequence;
    }

    public boolean setMer(int position, short value) {
        if (mers[position] == 0) {
            mers[position] = value;
            mersCount++;
            return true;
        }
        return false;
    }

    public int getMaxMers(int k) {
      return getSequence().getLength() - k + 1;
    }
    
    public int getMaxUniqMers() {
        return maxUniqMers;
    }

    public void setMaxUniqMers(int maxUniqMers) {
        this.maxUniqMers = maxUniqMers;
    }

    public short[] getMers() {
        return mers;
    }
    
    public void collectStatsAndResetMers(String sampleName, String TOOL_NAME) {
        samplesToStatsMap.put(sampleName, new KmerFilterStats(mers, mersCount, getMaxUniqMers()));
        
        mers = new short[mers.length];
        mersCount = 0;
    }

    public String getId() {
        return id;
    }

    public KmerFilterStats getKmerFilterStats(String sampleId) {
        return samplesToStatsMap.get(sampleId);
    }
    

    @Override
    public int compareTo(KmerFilter o) {
        return getId().compareTo(o.getId());
    }
    
    @Override
    public boolean equals(Object o) {
        return compareTo((KmerFilter) o) == 0;
    }

    @Override
    public int hashCode() {
        return getId().hashCode();
    }

    
}
