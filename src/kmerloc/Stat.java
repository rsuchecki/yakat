/*
 * Copyright 2020 rad.
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
package kmerloc;

/**
 *
 * @author rad
 */
public class Stat {
    private final String fastaID;
    private final int k;
    private final long count;

    public Stat(String fastaID, int k, long count) {
        this.fastaID = fastaID;
        this.k = k;
        this.count = count;
    }

    public String getFastaID() {
        return fastaID;
    }

    public int getK() {
        return k;
    }

    public long getCount() {
        return count;
    }
    
    
}
