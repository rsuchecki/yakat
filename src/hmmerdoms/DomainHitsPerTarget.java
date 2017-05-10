/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class DomainHitsPerTarget {

    private HashMap<String, ArrayList<DomainHit>> targetToHitsMap = new HashMap<>();
    private String SEP = "\t";

    public void addHit(DomainHit domainHit) {
        ArrayList<DomainHit> hits = targetToHitsMap.get(domainHit.getTargetId());
        if (hits == null) {
            hits = new ArrayList<>();
            targetToHitsMap.put(domainHit.getTargetId(), hits);
        }
        hits.add(domainHit);
    }

    public void printAll(boolean forward, int maxGap) {
        String[] sortedKeys = getSortedKeys();
        Integer prevEnd = null;
        String prevChr = null;
        int clustered = 0;
        Integer clusteredFrom = null;
        Integer clusteredTo = null;
        String clusteredChr = null;
        String clusteredStrand = null;
        Integer clusterElems = null;

        for (int i = 0; i < sortedKeys.length; i++) {
            ArrayList<DomainHit> hits = targetToHitsMap.get(sortedKeys[i]);
            if (forward) {
                Collections.sort(hits, new HitForwardComparator());
            } else {
                Collections.sort(hits, new HitReverseComparator());
            }
            for (DomainHit hit : hits) {
                String chr = hit.getTargetId().replaceAll("_.*", "");
                StringBuilder sb = new StringBuilder(chr);
                int frame = hit.getTargetFrame();
                int contigOffset = hit.getTargetOffset();
                int envStartLocal = -1;
                int envStart = -1;
                int envEndLocal = -1;
                int envEnd = -1;
                boolean print = false;
                String strand = "plus";
                if (forward && frame > 0) {
                    envStartLocal = hit.getEnvFrom() * 3 - (3-frame);
                    envStart = contigOffset + envStartLocal; 
                    envEndLocal = hit.getEnvTo() * 3 + (frame-1);
                    envEnd = contigOffset + envEndLocal;
                    print = true;
                } else if (!forward && frame < 0) {
                    int contigLen = hit.getTargetLength();
                    //TODO the adjustmant for frame might be off
                    envStartLocal = (contigLen - hit.getEnvTo() + 1) * 3 - (3+frame);
                    envEndLocal = (contigLen - hit.getEnvFrom() + 1) * 3 + (-frame-1);
                    envStart = contigOffset + envStartLocal; 
                    envEnd = contigOffset + envEndLocal;
                    print = true;
                    strand = "minus";
                }
                if (print) {
                    sb.append(SEP).append(envStart);
                    sb.append(SEP).append(envEnd);
                    sb.append(SEP).append(strand);
                    if (prevEnd != null && prevChr.equals(chr)) {
                        sb.append(SEP).append(envStart - prevEnd);
                    }
                    //printing domain info
                    System.err.println(sb.toString() + "\t" + hit.getTargetId()+ "\t"+envStartLocal+ "\t"+envEndLocal+ "\t"+hit.getEnvFrom()+ "\t"+hit.getEnvTo());
//                if(distanceToLast < 500 && clusteredStrand == null || clusteredStrand ==strand) {

                    if ((prevChr == null || prevChr.equals(chr)) && (prevEnd == null || envStart - prevEnd < maxGap)) {
                        clusteredFrom = clusteredFrom == null ? envStart : clusteredFrom > envStart ? envStart : clusteredFrom;
                        clusteredTo = clusteredTo == null ? envEnd : clusteredTo > envEnd ? envEnd : clusteredTo;
                        clusteredStrand = clusteredStrand == null ? strand : clusteredStrand.equals(strand) ? strand : "ERROR_CLUSTERED_ON_BOTH_STARNDS";
                        String tmp = "" + clusteredChr;
                        clusteredChr = clusteredChr == null ? chr : clusteredChr.equals(chr) ? chr : "ERROR_CLUSTERED_ACROSS_CHROMOSOMES";
                        if (clusteredChr.startsWith("ERROR")) {
                            int x = 0;
                        }
                        clusterElems = clusterElems == null ? 1 : ++clusterElems;
                    } else {
                        if (clusteredFrom != null) {
                            //printing group/cluster info
                            System.out.println(clusteredChr + SEP + clusteredFrom + SEP + clusteredTo + SEP + clusteredStrand + SEP + clusterElems);
                        }
                        clusteredFrom = null;
                        clusteredTo = null;
                        clusteredStrand = null;
                        clusteredChr = null;
                        clusterElems = null;
                    }
                    prevEnd = envEnd;
                    prevChr = chr;
                }
            }
        }
    }

    private String[] getSortedKeys() {
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
            return dh1.getEnvFrom() - dh2.getEnvFrom();
        }
    }

    private class HitReverseComparator implements Comparator<DomainHit> {

        @Override
        public int compare(DomainHit dh1, DomainHit dh2) {
            return dh2.getEnvFrom() - dh1.getEnvFrom();
        }
    }

}
