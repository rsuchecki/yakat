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

import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;
import shared.Sequence;

/**
 *
 * @author rad
 */
public class SeedSequence extends Sequence {
    ConcurrentHashMap<Integer, String> kToExtendedSequence;
    
    public SeedSequence(String id, String sequenceString) {
        super(id, sequenceString);        
        kToExtendedSequence = new ConcurrentHashMap<>();
    }
    
    public void setExtended(int k, String extendedSequenceString) {
        String previous = kToExtendedSequence.putIfAbsent(k, extendedSequenceString);
        if(previous != null) {
            Reporter.report("[ERROR]", "Unexpected value present at k="+k, this.getClass().getSimpleName());
        }
    }
    
    public String getExtended(int k) {
        return kToExtendedSequence.get(k);
    }
    
}
