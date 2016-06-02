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
package vsearchprocess;

/**
 *
 * @author rad
 */
public class Snp {

    private MsaSequence s1;
    private MsaSequence s2;
    private int position;

    /**
     *
     * @param s1
     * @param s2
     * @param position (zero indexed)
     */
    public Snp(MsaSequence s1, MsaSequence s2, int position) {
        this.s1 = s1;
        this.s2 = s2;
        this.position = position;
    }

    public MsaSequence getS1() {
        return s1;
    }

    public MsaSequence getS2() {
        return s2;
    }

    /**
     *
     * @return zero-indexed position
     */
    public int getPosition_0() {
        return position;
    }

    /**
     *
     * @return 1-indexed position
     */
    public int getPosition_1() {
        return position + 1;
    }

    public CharSequence getSnpString(int clusterNumber, boolean reverseLex, String DELIMITER, String suffix) {
        StringBuilder sb = new StringBuilder("Cluster_");
        sb.append(clusterNumber);
        String id1 = s1.getId();
        String id2 = s2.getId();
        int order = id1.compareTo(id2);
        order = reverseLex ? -order : order;
        if (order > 0) {
            sb.append(DELIMITER).append(id1).append(DELIMITER).append(s1.getLength());
            sb.append(DELIMITER).append(id2).append(DELIMITER).append(s2.getLength());
            sb.append(DELIMITER).append(getPosition_1()).append(DELIMITER);
            sb.append(s1.getSequenceString().charAt(getPosition_0())).append(DELIMITER).append(s2.getSequenceString().charAt(getPosition_0()));
        } else {
            sb.append(DELIMITER).append(id2).append(DELIMITER).append(s2.getLength());
            sb.append(DELIMITER).append(id1).append(DELIMITER).append(s1.getLength());
            sb.append(DELIMITER).append(getPosition_1()).append(DELIMITER);
            sb.append(s2.getSequenceString().charAt(getPosition_0())).append(DELIMITER).append(s1.getSequenceString().charAt(getPosition_0()));
        }
        sb.append(DELIMITER);
        if(suffix != null) {
            sb.append(suffix);
        }
        return sb;
    }
}
