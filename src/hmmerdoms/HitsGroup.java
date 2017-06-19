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

import java.text.NumberFormat;
import java.util.ArrayList;
import shared.Orf;
import shared.Sequence;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class HitsGroup {
    private ArrayList<DomainHit> domainHits;   
    private int from;
    private int to;
    private Strand strand;
    private enum Strand {
        plus, minus
    }
    
    
//    public HitsGroup(ArrayList<DomainHit> domainHits, String targetId) {
//        this.domainHits = domainHits;        
//        this.targetId = domainHits.get(0).getTargetIdCorrected();
//        
//        this.from = from;
//        this.to = to;
//        Sequence sequence = new Sequence(targetId+":"+from+"-"+to, fasta.getSequence(chr,clusteredFrom.longValue(),clusteredTo.longValue()));
//                            ArrayList<Orf> orfs = sequence.getOrfs();
//                            for (Orf orf : orfs) {
//        this.orfs = orfs;
//    }

    public HitsGroup(boolean forward) {
        this.strand = forward ? Strand.plus : Strand.minus;
    }

    
    
    public void addDomainHit(DomainHit domainHit) {
        int start = domainHit.getCorrectedStart();;
        int end = domainHit.getCorrectedEnd();
        if(domainHits == null) {
            domainHits = new ArrayList<>();
            from = start;
            to = end;
        } else {
            from = from > start ? start : from;
            to = to < end ? end : to;
        }
        domainHits.add(domainHit);
                
    }

    public ArrayList<DomainHit> getDomainHits() {
        return domainHits;
    }


    public String getTargetId() {
        return getDomainHits().get(0).getTargetIdCorrected();
    }

    public int getFrom() {
        return from;
    }

    public int getTo() {
        return to;
    }
    
    public int size() {
        return domainHits == null ? 0 : domainHits.size();
    }
    
    public boolean isEmpty() {
        return size() == 0;
    }

    public String getStrand() {
        return strand.toString();
    }
    
    
    public CharSequence getSummary(String sep) {
        StringBuilder sb = new StringBuilder(getTargetId());
        sb.append(sep).append(getFrom());
        sb.append(sep).append(getTo());
        sb.append(sep).append(NumberFormat.getInstance().format(getTo()-getFrom()+1));
        sb.append(sep).append(getStrand());
        sb.append(sep).append(size());
        sb.append(sep);
        for (DomainHit domainHit : domainHits) {
            sb.append(domainHit.getTargetFrame()).append(",");
        }
        return sb.substring(0, sb.length()-1);
    }
    
    
    
    
    
}
