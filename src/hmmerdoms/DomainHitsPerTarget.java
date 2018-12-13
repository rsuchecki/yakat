/*
 * Copyright 2017 Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au.
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
package hmmerdoms;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;
import shared.FastaIndexed;
import shared.Orf;
import shared.Sequence;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class DomainHitsPerTarget {

    private HashMap<String, ArrayList<DomainHit>> targetToHitsMap = new HashMap<>();
    private String SEP = "\t";

    /**
     * Hits are stored in a hashmap under a chromosome id (chr1A,chr1B,...,chr7D) except for chromosome "chrUn" 
     * for which the hits are stored separately for each contig to avoid domains being grouped due to spurious 
     * proximity on the artificial pseudomolecule
     * @param domainHit 
     */
    public void addHit(DomainHit domainHit) {
        String targetIdCorrected = domainHit.getTargetIdCorrected();        
        String targetId = targetIdCorrected.startsWith("chrUn") ? domainHit.getTargetId() : targetIdCorrected;
        ArrayList<DomainHit> hits = targetToHitsMap.get(targetId);
        if (hits == null) {
            hits = new ArrayList<>();
            targetToHitsMap.put(targetId, hits);
        }
        hits.add(domainHit);
    }

    public ArrayList<HitsGroup> processDomainsPerStrand(boolean forward, int maxGap, FastaIndexed fasta) {
        String[] sortedKeys = getSortedTargetIds();
        ArrayList<HitsGroup> hitsGroups = new ArrayList<>();
        for (int i = 0; i < sortedKeys.length; i++) {
//            System.err.println("Processing "+sortedKeys[i]);
            Integer prevEnd = null;
            ArrayList<DomainHit> hits = targetToHitsMap.get(sortedKeys[i]);
            //Strands processed separately
//            if (forward) {
                Collections.sort(hits, new HitForwardComparator());
//            } else {
//                Collections.sort(hits, new HitReverseComparator());
//            }
            HitsGroup hitsGroup = new HitsGroup(forward);
            for (DomainHit hit : hits) {
                int frame = hit.getTargetFrame();
                //Strands processed separately
                if ((forward && frame > 0) || (!forward && frame < 0)) {
                    int correctedStart = hit.getCorrectedStart();
                    int correctedEnd = hit.getCorrectedEnd();
                    

//                    System.err.println("dist = "+NumberFormat.getInstance().format( correctedStart - prevEnd));
                    //If this is the first record or within maximum gap distance to the last one                    
                    if (prevEnd != null &&  correctedStart - prevEnd > maxGap) {
//                        System.err.println("Adding "+hitsGroup.getDomainHits().size()+" as distance = "+( NumberFormat.getInstance().format(hit.getCorrectedStart() - prevEnd))+" at "+ hitsGroup.getTargetId());
                        hitsGroups.add(hitsGroup);
                        hitsGroup = new HitsGroup(forward);
                    }
                    hitsGroup.addDomainHit(hit);
//                    System.err.println(sb.toString() + "\t" + hit.getTargetId() + "\t" + startLocal + "\t" + endLocal + "\t" + hit.getAlnFrom() + "\t" + hit.getAlnTo());
                    prevEnd = hit.getCorrectedEnd();
                }
            }
            if (!hitsGroup.isEmpty()) {
                hitsGroups.add(hitsGroup);
//                    System.err.println("Adding " + hitsGroup.getDomainHits().size() + " as last one " + hitsGroup.getTargetId());
            }
        }
        return hitsGroups;
    }

    private String[] getSortedTargetIds() {
        Set<String> keySet = targetToHitsMap.keySet();
        String keys[] = new String[keySet.size()];
        keys = keySet.toArray(keys);
        Arrays.sort(keys, new KeyComparator());
        return keys;
    }

    private class KeyComparator implements Comparator<String> {

        @Override
        public int compare(String o1, String o2) {
            String[] split1 = o1.split("_");
            String[] split2 = o2.split("_");
            int comp = split1[0].compareTo(split2[0]);
            if (comp == 0) {
                Integer offset1 = Integer.parseInt(split1[1]);
                Integer offset2 = Integer.parseInt(split2[1]);
                return offset1.compareTo(offset2);
            } else {
                return comp;
            }
        }

    }

    private class HitForwardComparator implements Comparator<DomainHit> {

        @Override
        public int compare(DomainHit dh1, DomainHit dh2) {
            return dh1.getCorrectedStart()- dh2.getCorrectedStart();
        }
    }

    private class HitReverseComparator implements Comparator<DomainHit> {

        @Override
        public int compare(DomainHit dh1, DomainHit dh2) {
            return dh2.getCorrectedStart() - dh1.getCorrectedStart();
        }
    }

}
