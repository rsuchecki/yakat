/*
 * Copyright 2020 rad.
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
package kmermatch;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import kmermatch.Kmer;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class KmerSetsMap {

    private final ConcurrentHashMap<Integer, ConcurrentSkipListSet<Kmer>> kmersSetsMap;
    private final ArrayList<Integer> kSizes;
    private final String TOOL_NAME;
    
    public KmerSetsMap(String TOOL_NAME) {
        this.kmersSetsMap = new ConcurrentHashMap<>();
        this.kSizes = new ArrayList<>();
        this.TOOL_NAME = TOOL_NAME;
    }

    public ConcurrentHashMap<Integer, ConcurrentSkipListSet<Kmer>> getKmerSetsMap() {
        return kmersSetsMap;
    }

    public ConcurrentSkipListSet<Kmer> getKmerSet(int k) {
        ConcurrentSkipListSet<Kmer> set = kmersSetsMap.get(k);
        if (set == null) {
            set = new ConcurrentSkipListSet<>();
            ConcurrentSkipListSet previous = kmersSetsMap.putIfAbsent(k, set);
            if (previous == null) {
                addKValue(k);
            } else { //another thread just beat us to putting this one in
                return previous;
            }
        }
        return set;
    }
    
        private synchronized void addKValue(int k) {
        if (!kSizes.contains(k)) { //might have been added by input reader, 
            Reporter.report("[INFO]", "Initiating set for previously unseen k-mers, k=" + k, TOOL_NAME);
            kSizes.add(k);
        }
    }

}
