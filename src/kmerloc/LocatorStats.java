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
package kmerloc;

import java.util.Collection;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author rad
 */
public class LocatorStats {
    private final ConcurrentHashMap<String, ConcurrentHashMap<Integer, Integer>> statsMap;

    public LocatorStats() {
        statsMap = new ConcurrentHashMap<>();
    }
   
    public synchronized void addStat(String fastaId, int k, int count) {
        ConcurrentHashMap<Integer, Integer> statMap = statsMap.get(fastaId);
        if(statMap == null) {
            statMap = new ConcurrentHashMap<>();
            statsMap.put(fastaId, statMap);            
        }       
        Integer storedCount = statMap.get(k);
        if(storedCount == null) {
            statMap.put(k, count);
        } else {
            statMap.put(k, count+storedCount);
        }
        
    }
    
    public int getStat(String fastaId, int k) {
        return statsMap.get(fastaId).get(k);
    }
    
    public int getTotalPlaced() {        
        return statsMap.values()
                .stream()
                .map(Map::values)
                .flatMap(Collection::stream)
                .reduce(0, Integer::sum);
    }
    
//    public int getPlacedPerK(int k) {
//        
//    }
       
    public int getPlacedPerId(String id) {        
        return statsMap.get(id).values().stream().reduce(0, Integer::sum);

        
    }
    
   
}
