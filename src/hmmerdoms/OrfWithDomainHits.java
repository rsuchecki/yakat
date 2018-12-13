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

import java.util.ArrayList;
import shared.Orf;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class OrfWithDomainHits implements Comparable<OrfWithDomainHits>{
    private Orf orf;
    private ArrayList<DomainHit> domainHits;

    public OrfWithDomainHits(Orf orf, ArrayList<DomainHit> domainHits) {
        this.orf = orf;
        this.domainHits = domainHits;
    }

    public Orf getOrf() {
        return orf;
    }

    public ArrayList<DomainHit> getDomainHits() {
        return domainHits;
    }

    @Override
    public int compareTo(OrfWithDomainHits o) {
        return (int) (orf.getFrom() - o.getOrf().getFrom());
    }
    

}
