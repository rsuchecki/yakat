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
import java.util.Collections;
import java.util.regex.Pattern;
import shared.FastaIndexed;
import shared.Orf;
import shared.Sequence;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class HitsGroup implements Comparable<HitsGroup>{

    private ArrayList<DomainHit> domainHits;
    private int from;
    private int to;
    private Long orfPredictionRangeFrom;
    private Long orfPredictionToPosition;
    private Strand strand;

    ArrayList<OrfWithDomainHits> orfsWithDomainsHits;

    @Override
    public int compareTo(HitsGroup o) {
        int lexId = getTargetId().compareTo(o.getTargetId());
        if(lexId == 0) {
            return getFrom() - o.getFrom();
        } else {
            return lexId;
        }
    }

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
        if (domainHits == null) {
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
    
    public String getStrandSymbol() {
        return strand == Strand.plus ? "+" : "-";
    }

    public CharSequence getSummary(String sep) {
        StringBuilder sb = new StringBuilder(getTargetId());
        sb.append(sep).append(getFrom());
        sb.append(sep).append(getTo());
        sb.append(sep).append(NumberFormat.getInstance().format(getTo() - getFrom() + 1));
        sb.append(sep).append(getStrand());
        sb.append(sep).append(size());
        sb.append(sep);
        for (DomainHit domainHit : domainHits) {
            sb.append(domainHit.getTargetFrame()).append(",");
        }
        return sb.substring(0, sb.length() - 1);
    }

    public int getNumAssigned() {
        int assigned = 0;
        for (DomainHit domainHit : domainHits) {
            if (domainHit.isAssigned()) {
                assigned++;
            }
        }
        return assigned;
    }

    private void identifyOrfsOverlappingDomains(FastaIndexed fastaIndexed, int upstreamDownstreamBases, boolean allowRecursiveExtending) {
        orfsWithDomainsHits = new ArrayList<>();
        orfPredictionRangeFrom = new Long(getFrom()) - upstreamDownstreamBases;
//        int mod = (int) (orfPredictionFromPosition % 3);
//        switch(mod) {
//            case 0: orfPredictionFromPosition -= 2; break;
//            case 1: break;
//            case 2: orfPredictionFromPosition -= 1; break;
//        }
//        orfPredictionFromPosition -= mod+1;
        orfPredictionToPosition = new Long(getTo()) + upstreamDownstreamBases;
        Sequence s = new Sequence(getTargetId(), fastaIndexed.getSequence(getTargetId(), orfPredictionRangeFrom, orfPredictionToPosition));
                Pattern startStopCodonsForward = Pattern.compile("ATG|(T(AA|AG|GA))", Pattern.CASE_INSENSITIVE);
        Pattern startStopCodonsReverse = Pattern.compile("CAT|((TT|TC|CT)A)", Pattern.CASE_INSENSITIVE);
        Pattern stopCodonsForward = Pattern.compile("T(AA|AG|GA)", Pattern.CASE_INSENSITIVE);
        Pattern stopCodonsReverse = Pattern.compile("(TT|TC|CT)A", Pattern.CASE_INSENSITIVE);
        
        ArrayList<Orf> orfs = s.getOrfs(0, startStopCodonsForward, startStopCodonsReverse, stopCodonsForward, stopCodonsReverse);
        Collections.sort(orfs);
        boolean extendFurther = false;
        for (Orf orf : orfs) {
            ArrayList<DomainHit> domainsInOrf = new ArrayList<>();
            long orfFrom = orf.getFrom() + orfPredictionRangeFrom - 1;
            long orfTo = orf.getTo() + orfPredictionRangeFrom - 1;
//                System.err.println("ORF:" + orf.getParenId() + "\t" + orf.getFrom() + "\t" + orf.getTo() + "\t" + orf.getLength() + "\t" + orfFrom + "\t" + orfTo + "\t" + orf.getFrame() + "\t" + orf.hasStopCodon());
            int orfFrame = orf.getFrame();
            for (DomainHit hit : getDomainHits()) {
                int hitFrame = hit.getTargetFrame();
                if ((orfFrame > 0 && hitFrame > 0) || (orfFrame < 0 && hitFrame < 0)) {
                    int hitFrom = hit.getCorrectedStart();
                    int hitTo = hit.getCorrectedEnd();
                    //IDENTIFY DOMAINS CONTAINED IN OR OVERLAPPING WITH THE ORF
                    if ((orfFrom <= hitFrom && orfTo >= hitTo) || (orfFrom >= hitFrom && orfFrom <= hitTo) || (orfTo >= hitFrom && orfTo <= hitTo)) {
//                            domainsInOrf++;
                        if (!hit.isAssigned()) {
                            domainsInOrf.add(hit);
                            hit.setAssigned(true);
                        }
//                        } else if((orfFrom >= hitFrom && orfFrom <= hitTo) || (orfTo >= hitFrom && orfTo <= hitTo)) {
//                            domainsOverlapOrf.add(hit);
//                        } else {                            
//                            domainsOutOrf.add(hit);
                    }
                }
            }
            //IF a domain-containing ORF is unterminated, try extending             
            if (allowRecursiveExtending && domainsInOrf.size() > 0) {
                if (!orf.hasStopCodon()) {
                    extendFurther = true;
                    break;
                }
                orfsWithDomainsHits.add(new OrfWithDomainHits(orf, domainsInOrf));
//                System.err.println("ORF:    " + orf.getParenId() + "\t" + orfFrom + "\t" + orfTo + "\t" + orf.getFrame() + "\t" + (orf.hasStopCodon() ? "" : "non-stop"));
//                for (int i = 0; i < domainsInOrf.size(); i++) {
//                    DomainHit dh = domainsInOrf.get(i);
////                    System.err.println("DOM_in :" + dh.getTargetIdCorrected() + "\t" + dh.getCorrectedStart() + "\t" + dh.getCorrectedEnd() + "\t" + dh.getTargetFrame() + "\t<=== " + (i + 1));
//                }
//                System.err.println(s.getSequenceString().substring((int) orf.getFrom() - 1, Math.min(s.getLength(), (int) orf.getTo())));
            }
        }
        if (extendFurther) {
            for (DomainHit hit : getDomainHits()) {
                hit.setAssigned(false); //resret previous assignment to ORFs 
            }
//            System.err.println("Further extending +- "+(upstreamDownstreamBases*2));
            identifyOrfsOverlappingDomains(fastaIndexed, upstreamDownstreamBases*2, allowRecursiveExtending);
        }

    }

    public ArrayList<OrfWithDomainHits> getOrfsWithDomainsHits(FastaIndexed fastaIndexed, int upstreamDownstreamBases, boolean allowRecursiveExtending) {
        if (orfsWithDomainsHits == null) {
            identifyOrfsOverlappingDomains(fastaIndexed, upstreamDownstreamBases, allowRecursiveExtending);
        }
        return orfsWithDomainsHits;
    }

    public Long getOrfPredictionFromPosition() {
        return orfPredictionRangeFrom;
    }

    public Long getOrfPredictionToPosition() {
        return orfPredictionToPosition;
    }

    public CharSequence getGff3() {
        return "";
    }

}
