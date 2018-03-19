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
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import shared.Reporter;

/**
 * A sequence with a map kValue -> extension
 *
 * @author rad
 */
public class SeedSequence {

    private final String id;
    private final String sequenceString;
//    private final ConcurrentHashMap<Integer, String> kToExtendedSequence;
    private final ConcurrentHashMap<Integer, SeedExtensionsPair> kToExtensionsPair;

//    private final ConcurrentHashMap<Integer, Integer> kToExtensionCount;
//    private final int[] extensionCounts = new int[101];
//    private ConcurrentHashMap<Integer, PairMer> kToHeadMer;
//    private ConcurrentHashMap<Integer, PairMer> kToTailMer;
    public SeedSequence(String id, String sequenceString) {
        this.id = id;
//        kToExtendedSequence = new ConcurrentHashMap<>();
//        kToExtendedSequence.put(0, sequenceString);
        this.sequenceString = sequenceString;

        kToExtensionsPair = new ConcurrentHashMap<>();
        kToExtensionsPair.put(0, new SeedExtensionsPair());
//        kToHeadMer = new ConcurrentHashMap<>();
//        kToTailMer = new ConcurrentHashMap<>();
////        kToExtensionCount = new ConcurrentHashMap<>();
    }

    public SeedSequence() {
        this.id = null;
        this.sequenceString = null;
        this.kToExtensionsPair = null;
    }

    public void setLeftExtension(int k, String leftExtension) {
        SeedExtensionsPair pair = kToExtensionsPair.get(k);
        if (pair == null) {
            pair = new SeedExtensionsPair();
            kToExtensionsPair.put(k, pair);
        }
        pair.setExtensionLeft(leftExtension);
    }

    public void setRightExtension(int k, String rightExtension) {
        SeedExtensionsPair pair = kToExtensionsPair.get(k);
        if (pair == null) {
            pair = new SeedExtensionsPair();
            kToExtensionsPair.put(k, pair);
        }
        pair.setExtensionRight(rightExtension);
    }

//    public void setExtensions(int k, String leftExtension, String rightExtension) {
//        SeedExtensionsPair pair = kToExtensionsPair.get(k);
//        pair.setExtensionLeft(leftExtension);
//        pair.setExtensionRight(rightExtension);
//    }
//    public synchronized void setExtended(int k, String extendedSequenceString) { //, int extensions) {
//        if (extendedSequenceString.length() > getSequenceString().length()) {
//            String previous = kToExtendedSequence.putIfAbsent(k, extendedSequenceString);
//            if (previous != null) {
//                if (previous.length() < extendedSequenceString.length()) {
//                    kToExtendedSequence.put(k, extendedSequenceString);
////                    extensionCounts[k]++;
////                    System.err.println("at k=" + k + " current count = " + extensionCounts[k] +" updated ");
////                    System.err.println(extendedSequenceString);
//                }
////            Reporter.report("[ERROR]", "Unexpected value present at k=" + k, this.getClass().getSimpleName());
////            if (!previous.equals(extendedSequenceString)) {
////                System.err.println(previous + "<-pre");
////                System.err.println(extendedSequenceString + "<-new");
////            }
//            } else {
////                extensionCounts[k]++;
////                System.err.println("at k=" + k + " current count = " + extensionCounts[k] +" new ");
////                System.err.println(extendedSequenceString);
//            }
////        } else {
////            System.err.println("trying to store short extension - investigate "+this.getClass().getName());
////            
//        }
////        Integer previousCount = kToExtensionCount.putIfAbsent(k, extensions);
////        if(previousCount != null) {
////            kToExtensionCount.put(k, previousCount+extensions);
////        }
//    }
//    /**
//     * Get extended sequence for a given value of k
//     *
//     * @param k
//     * @return extended sequence or null
//     */
//    public String getExtended(int k) {
//        return kToExtendedSequence.get(k);
//    }
//    /**
//     * Get extended sequence for aSeedSequence given value of k, or the original sequence
//     *
//     * @param k
//     * @return extended sequence or null
//     */
//    public String getExtendedOrOriginal(int k) {
//        return kToExtendedSequence.containsKey(k) ? kToExtendedSequence.get(k) : kToExtendedSequence.get(0);
//    }
//    /**
//     *
//     * @return longest extended String or input if longest not longer than input original String
//     */
//    public Map.Entry<Integer, String> getLongestExtended() {
//        Map.Entry<Integer, String> longest = null;
//        Map.Entry<Integer, String> input = null;
//
////        int[] lens = new int[100];
//        for (Map.Entry<Integer, String> entry : kToExtendedSequence.entrySet()) {
//            if (entry.getKey() == 0) {
//                input = entry;
//            }
//            String value = entry.getValue();
//            longest = (longest == null || (longest.getValue().length() < value.length())) ? entry : longest;
////            lens[entry.getKey()] = entry.getValue().length();
//        }
////        System.err.println(id);
////        for (int i = 0; i < extensionCounts.length; i++) {
////            if (extensionCounts[i] > 0) {
////                System.err.println("[" + i + "] = " + extensionCounts[i]);
////            }
////        }
//
//        return longest.getValue().length() == getSequenceString().length() ? input : longest;
//
//    }
    /**
     * Iterate through per-k extensions and output the longest one. Iteration in no particular order so if multiple
     * values of k point to an equally long extension,
     *
     * @return
     */
    public Map.Entry<Integer, SeedExtensionsPair> getLongestExtensionLeft(String TOOL_NAME) {
        Map.Entry<Integer, SeedExtensionsPair> longest = null;
        for (Map.Entry<Integer, SeedExtensionsPair> entry : kToExtensionsPair.entrySet()) {
            String value = entry.getValue().getExtensionLeft();
            longest = (longest == null || (longest.getValue().getExtensionLeft().length() < value.length())) ? entry : longest;
            if (longest.getValue().getExtensionLeft().length() == value.length() && !longest.getValue().getExtensionLeft().equals(value)) {
                String mesg = "Same length but different seed extensions, "
                    + longest.getValue().getExtensionLeft()+ "@k=" + longest.getKey() + " " + value + "@k="
                    + entry.getKey() + ", L, seed=" + getId();
                longest = longest.getKey() < entry.getKey() ? entry : longest; //store higher k extension
                mesg += ", storing extension @k=" + longest.getKey();
                Reporter.report("[WARNING]", mesg, TOOL_NAME);
            }
        }
        return longest;
    }

    public Map.Entry<Integer, SeedExtensionsPair> getLongestExtensionRight(String TOOL_NAME) {
        Map.Entry<Integer, SeedExtensionsPair> longest = null;
        for (Map.Entry<Integer, SeedExtensionsPair> entry : kToExtensionsPair.entrySet()) {
            String value = entry.getValue().getExtensionRight();
            longest = (longest == null || (longest.getValue().getExtensionRight().length() < value.length())) ? entry : longest;
            if (longest.getValue().getExtensionRight().length() == value.length() && !longest.getValue().getExtensionRight().equals(value)) {
                String mesg = "Same length but different seed extensions, "
                    + longest.getValue().getExtensionRight() + "@k=" + longest.getKey() + " " + value + "@k="
                    + entry.getKey() + ", R, seed=" + getId();
                longest = longest.getKey() < entry.getKey() ? entry : longest; //store higher k extension
                mesg += ", storing extension @k=" + longest.getKey();
                Reporter.report("[WARNING]", mesg, TOOL_NAME);
            }
        }
        return longest;
    }

    public String getId() {
        return id;
    }

    public String getSequenceString() {
        return sequenceString;
//        return kToExtendedSequence.get(0);
    }

//    public void generateEndPairMers(int k) {
//            PairMer headMer = PairMerGenerator.generatePairMer(getSequenceString().substring(0, k), false, k - 1);
//            int len = getSequenceString().length();
//            PairMer tailMer = PairMerGenerator.generatePairMer(getSequenceString().substring(len - k, len), true, k - 1);
//            kToHeadMer.putIfAbsent(k, headMer);
//            kToTailMer.putIfAbsent(k, tailMer);
//    }
}
