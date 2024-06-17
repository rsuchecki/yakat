/*
 * Copyright 2017 rad.
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

import java.util.ArrayList;

/**
 *
 * @author rad
 */
public class KmerFilterStats {

//    private final short[] mers;
    private final double medianFrequency;
    private final int mersCount;
    private final int maxUniqMers;

    public KmerFilterStats(short[] mers, int mersCount, int maxUniqMers) {
//        this.mers = mers;
        this.medianFrequency = shared.CommonMaths.getMedian(getNonZeroKmerFreqs(mers));
        this.mersCount = mersCount;
        this.maxUniqMers = maxUniqMers;
    }

    private ArrayList<Short> getNonZeroKmerFreqs(short[] mers) {
        ArrayList<Short> nonZeroValues = new ArrayList<>(mers.length);
        for (Short value : mers) {
            if (value > 0) {
                nonZeroValues.add(value);
            }
        }
        return nonZeroValues;
    }

//    private double getMax(ArrayList<Short> nonZeroValues) {
//        if (nonZeroValues.isEmpty()) {
//            return 0;
//        } else {
//            short max = 0;
//            for (Short nonZeroValue : nonZeroValues) {
//                max = nonZeroValue > max ? nonZeroValue : max;
//            }
//            return max;
//        }
//    }

    public int getMersCount() {
        return mersCount;
    }

    public int getMaxUniqMers() {
        return maxUniqMers;
    }
    
//    public int getMaxMers() {
//      return mers.length;
//    }

  public double getMedianFrequency() {
    return medianFrequency;
  }
    
    
    
//    public double getMedianFrequency() {
//        return shared.CommonMaths.getMedian(getNonZeroKmerFreqs(mers));
//    }

//    public int getCoverage() {
//        return getNonZeroKmerFreqs(mers).size();
//    }

    public double getCoverageRatioOverMaxUniqmers() {
//        int coverage = getCoverage();
//        int maxMers1 = getMaxMers();
//        return (double) getCoverage() / (double) getMaxUniqMers();
        return (double) getMersCount() / (double) getMaxUniqMers();
    }

//    public short[] getMers() {
//        return mers;
//    }
    
    
}
