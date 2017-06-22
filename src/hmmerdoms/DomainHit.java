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

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class DomainHit {

    private String targetId; //1
    private int targetLength; //3
    private String queryId; //4
    private int queryLength; //6
    private double iEvalue; //13
    private double score; //14
    private double bias; //15
    private int profileFrom; //16
    private int profileTo; //17
    private int alnFrom; //18
    private int alnTo; //19
    private int envFrom; //20
    private int envTo; //21
    private double accPosterior;  //22  The mean posterior probability of aligned residues in the MEA alignment; a measure of how reliable  the  overall  alignment  is  (from  0  to  1,  with  1.00  indicating  a  completely  reliable  alignment according to the model).
    private String targetDescription; //23

    //custom
    private int targetFrame;
    private boolean assigned;


    public DomainHit(String[] domtblRecords) {
        for (int i = 0; i < domtblRecords.length; i++) {
            int pos = i + 1;
            switch (pos) {
                case 1:
                    targetId = domtblRecords[i];
                    break;
                case 3:
                    targetLength = Integer.parseInt(domtblRecords[i]);
                    break;
                case 4:
                    queryId = domtblRecords[i];
                    break;
                case 6:
                    queryLength = Integer.parseInt(domtblRecords[i]);
                    break;
                case 13:
                    iEvalue = Double.parseDouble(domtblRecords[i]);
                    break;
                case 14:
                    score = Double.parseDouble(domtblRecords[i]);
                    break;
                case 15:
                    bias = Double.parseDouble(domtblRecords[i]);
                    break;
                case 16:
                    profileFrom = Integer.parseInt(domtblRecords[i]);
                    break;
                case 17:
                    profileTo = Integer.parseInt(domtblRecords[i]);
                    break;
                case 18:
                    alnFrom = Integer.parseInt(domtblRecords[i]);
                    break;
                case 19:
                    alnTo = Integer.parseInt(domtblRecords[i]);
                    break;
                case 20:
                    envFrom = Integer.parseInt(domtblRecords[i]);
                    break;
                case 21:
                    envTo = Integer.parseInt(domtblRecords[i]);
                    break;
                case 22:
                    accPosterior = Double.parseDouble(domtblRecords[i]);
                    break;
                case 23:
                    targetDescription = domtblRecords[i];
                    // Assuming this can be picked up from the description line for the target,
                    // BBMap translate6frames.sh appends fr{1,2,...,6} after a whitespace
                    int frame = Integer.parseInt(targetDescription.replaceFirst("fr", ""));
                    switch (frame) {
                        case 4:
                            targetFrame = -1;
                            break;
                        case 5:
                            targetFrame = -2;
                            break;
                        case 6:
                            targetFrame = -3;
                            break;
                        default:
                            targetFrame = frame;
                            break;
                    }
                    break;
            }

        }
    }

    public String getTargetId() {
        return targetId;
    }

    public int getTargetLength() {
        return targetLength;
    }

    public String getQueryId() {
        return queryId;
    }

    public int getQueryLength() {
        return queryLength;
    }

    public double getiEvalue() {
        return iEvalue;
    }

    public double getScore() {
        return score;
    }

    public double getBias() {
        return bias;
    }

    public int getProfileFrom() {
        return profileFrom;
    }

    public int getProfileTo() {
        return profileTo;
    }

    public int getAlnFrom() {
        return alnFrom;
    }

    public int getAlnTo() {
        return alnTo;
    }

    public int getEnvFrom() {
        return envFrom;
    }

    public int getEnvTo() {
        return envTo;
    }

    public double getAccPosterior() {
        return accPosterior;
    }

    public String getTargetDescription() {
        return targetDescription;
    }

    public int getTargetFrame() {
        return targetFrame;
    }

    /**
     * Assuming this can be picked up from the identifier, expected format:
     * >chromosome_contigOffset_Ns_length
     *
     * @return
     */
    public int getTargetOffset() {
        String[] toks = getTargetId().split("_");
        return Integer.parseInt(toks[1]);
    }

    public int getStartLocal() {
        int frame = getTargetFrame();
        if (frame > 0) {
            return getAlnFrom() * 3 - (3 - frame);
        } else {
            //TODO the adjustmant for frame might be off
            return (getTargetLength() - getAlnTo() + 1) * 3 - (3 + frame);
        }
    }

    public int getEndLocal() {
        int frame = getTargetFrame();
        if (frame > 0) {
            return getAlnTo() * 3 + (frame - 1);
        } else {
            return (getTargetLength() - getAlnFrom() + 1) * 3 + (-frame - 1);
        }
    }

    public Integer getCorrectedStart() {
        return getTargetOffset() + getStartLocal();
    }

    public Integer getCorrectedEnd() {
        return getTargetOffset() + getEndLocal();
    }

    public String getTargetIdCorrected() {
        return getTargetId().replaceAll("_.*", "");
    }

    public boolean isAssigned() {
        return assigned;
    }

    public void setAssigned(boolean assigned) {
        this.assigned = assigned;
    }
    
    
}

//From HMMER MANUAL, 1-indexed
//(1) target name: The name of the target sequence or profile.
//(2) target accession: Accession of the target sequence or profile, or ’-’ if none is available.
//(3) tlen: Length of the target sequence or profile, in residues.  This (together with the query length) is
//useful for interpreting where the domain coordinates (in subsequent columns) lie in the sequence.
//(4) query name: Name of the query sequence or profile.
//(5) accession: Accession of the target sequence or profile, or ’-’ if none is available.
//(6) qlen: Length of the query sequence or profile, in residues.
//(7) E-value: E-value of the overall sequence/profile comparison (including all domains).
//(8) score: Bit score of the overall sequence/profile comparison (including all domains),  inclusive of a
//null2 bias composition correction to the score.
//(9) bias: The biased composition score correction that was applied to the bit score.
//(10) #: This domain’s number (1..ndom).
//(11) of: The total number of domains reported in the sequence, ndom.
//(12) c-Evalue: The “conditional E-value”, a permissive measure of how reliable this particular domain
//may be.  The conditional E-value is calculated on a smaller search space than the independent E-
//value.  The conditional E-value uses the number of targets that pass the reporting thresholds.  The
//null hypothesis test posed by the conditional E-value is as follows.   Suppose that we believe that
//there is already sufficient evidence (from other domains) to identify the set of reported sequences as
//homologs of our query; now, how many additional domains would we expect to find with at least this
//particular domain’s bit score, if the rest of those reported sequences were random nonhomologous
//sequence (i.e. outside the other domain(s) that were sufficient to identified them as homologs in the first place)?
//(13) i-Evalue: The “independent E-value”, the E-value that the sequence/profile comparison would have
//received if this were the only domain envelope found in it, excluding any others.  This is a stringent
//measure  of  how  reliable  this  particular  domain  may  be.   The  independent  E-value  uses  the  total
//number of targets in the target database.
//(14) score: The bit score for this domain.
//(15) bias: The biased composition (null2) score correction that was applied to the domain bit score.
//(16) from (hmm coord): The start of the MEA alignment of this domain with respect to the profile, num-
//bered 1..N for a profile of N consensus positions.
//(17) to (hmm coord): The end of the MEA alignment of this domain with respect to the profile, numbered
//1..N for a profile of N consensus positions.
//(18) from (ali coord): The start of the MEA alignment of this domain with respect to the sequence,
//numbered 1..L for a sequence of L residues.
//(19) to (ali coord): The end of the MEA alignment of this domain with respect to the sequence, num-
//bered 1..L for a sequence of L residues.
//(20) from (env coord): The start of the domain envelope on the sequence, numbered 1..L for a se-
//quence of L residues.  The envelope defines a subsequence for which their is substantial probability
//mass supporting a homologous domain, whether or not a single discrete alignment can be identified.
//The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for
//weakly scoring domains.
//(21) to (env coord): The end of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
//(22) acc: The mean posterior probability of aligned residues in the MEA alignment; a measure of how
//reliable  the  overall  alignment  is  (from  0  to  1,  with  1.00  indicating  a  completely  reliable  alignment
//according to the model).
//(23) description of target: The remainder of the line is the target’s description line, as free text.
