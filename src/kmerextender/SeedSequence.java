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

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * A sequence with a map kValue -> extension
 *
 * @author rad
 */
public class SeedSequence {

    private final String id;
    private final ConcurrentHashMap<Integer, String> kToExtendedSequence;
//    private final ConcurrentHashMap<Integer, Integer> kToExtensionCount;
    private int[] extensionCounts = new int[101];

    public SeedSequence(String id, String sequenceString) {
        this.id = id;
        kToExtendedSequence = new ConcurrentHashMap<>();
        kToExtendedSequence.put(0, sequenceString);
//        kToExtensionCount = new ConcurrentHashMap<>();
    }

    public synchronized void setExtended(int k, String extendedSequenceString) { //, int extensions) {
        if (extendedSequenceString.length() > getSequenceString().length()) {
            String previous = kToExtendedSequence.putIfAbsent(k, extendedSequenceString);
            if (previous != null) {
                if (previous.length() < extendedSequenceString.length()) {
                    kToExtendedSequence.put(k, extendedSequenceString);
                    extensionCounts[k]++;
                    System.err.println("at k=" + k + " current count = " + extensionCounts[k]);
//                    System.err.println(extendedSequenceString);
                }
//            Reporter.report("[ERROR]", "Unexpected value present at k=" + k, this.getClass().getSimpleName());
//            if (!previous.equals(extendedSequenceString)) {
//                System.err.println(previous + "<-pre");
//                System.err.println(extendedSequenceString + "<-new");
//            }
            } else {
                extensionCounts[k]++;
                System.err.println("at k=" + k + " current count = " + extensionCounts[k]);
//                System.err.println(extendedSequenceString);
            }
        }
//        Integer previousCount = kToExtensionCount.putIfAbsent(k, extensions);
//        if(previousCount != null) {
//            kToExtensionCount.put(k, previousCount+extensions);
//        }
    }

    /**
     * Get extended sequence for a given value of k
     *
     * @param k
     * @return extended sequence or null
     */
    public String getExtended(int k) {
        return kToExtendedSequence.get(k);
    }

    /**
     * Get extended sequence for a given value of k, or the original sequence
     *
     * @param k
     * @return extended sequence or null
     */
    public String getExtendedOrOriginal(int k) {
        return kToExtendedSequence.containsKey(k) ? kToExtendedSequence.get(k) : kToExtendedSequence.get(0);
    }

    /**
     *
     * @return longest extended String or input if longest not longer than input
     * original String
     */
    public Map.Entry<Integer, String> getLongestExtended() {
        Map.Entry<Integer, String> longest = null;
        Map.Entry<Integer, String> input = null;

//        int[] lens = new int[100];
        for (Map.Entry<Integer, String> entry : kToExtendedSequence.entrySet()) {
            if (entry.getKey() == 0) {
                input = entry;
            }
            String value = entry.getValue();
            longest = (longest == null || (longest.getValue().length() < value.length())) ? entry : longest;
//            lens[entry.getKey()] = entry.getValue().length();
        }
        System.err.println(id);
        for (int i = 0; i < extensionCounts.length; i++) {
            if (extensionCounts[i] > 0) {
                System.err.println("[" + i + "] = " + extensionCounts[i]);
            }
        }

        return longest.getValue().length() == getSequenceString().length() ? input : longest;

    }

    public String getId() {
        return id;
    }

    public String getSequenceString() {
        return kToExtendedSequence.get(0);
    }
}
