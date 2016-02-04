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
package kmerextender;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 *
 * @author rad
 */
public class PairMerMaps {

    private final ConcurrentHashMap<Integer, PairMersMap> kSizeToPairMersMap;
    private final ArrayList<Integer> kSizes;

    public PairMerMaps(ArrayList<Integer> kSizes) {
        this.kSizes = kSizes;
        kSizeToPairMersMap = new ConcurrentHashMap<>(kSizes.size());
        for (Integer k : kSizes) {
            if (k != 0) {
                kSizeToPairMersMap.put(k, new PairMersMap(k));
            }
        }
    }

    public PairMersMap getPairMersMap(int k) {
        PairMersMap map = kSizeToPairMersMap.get(k);
        if(map == null) {
            map = new PairMersMap(k);
            kSizeToPairMersMap.put(k, map);
        }
        return map;
    }

    /**
     * Add a new PairMerMap for the given k
     *
     * @param k
     * @return false if map for given k already in
     */
    public boolean addK(int k) {
        if (kSizeToPairMersMap.putIfAbsent(k, new PairMersMap(k)) == null) {
            kSizes.add(k);
            return true;
        }
        return false;
    }

    /**
     * Add new PairMerMaps for the given k sizes
     *
     * @param kSizes
     * @return false if map for given k already in
     */
    public boolean addK(ArrayList<Integer> kSizes) {
        boolean kAlreadyIn = false;
        for (Integer k : kSizes) {
            if (kSizeToPairMersMap.putIfAbsent(k, new PairMersMap(k)) == null) {
                kSizes.add(k);
            } else {
                kAlreadyIn = true;
            }
        }
        return !kAlreadyIn;
    }

    public ArrayList<Integer> getkSizes() {
        return kSizes;
    }

    public int size() {
        if (kSizes.size() != kSizeToPairMersMap.size()) {
            Reporter.report("[BUG]", "Internal structures sizes mismtach", this.getClass().getCanonicalName());
        }
        return kSizes.size();
    }

}
