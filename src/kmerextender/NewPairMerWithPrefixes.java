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
package kmerextender;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class NewPairMerWithPrefixes {
    private int prefix;
    private int postfix;
    private PairMer pairMer;

    public NewPairMerWithPrefixes(int prefix, int postfix, PairMer pairMer) {
        this.prefix = prefix;
        this.postfix = postfix;
        this.pairMer = pairMer;
    }

    public int getPrefix() {
        return prefix;
    }

    public int getPostfix() {
        return postfix;
    }

    public PairMer getPairMer() {
        return pairMer;
    }
    
    
}
