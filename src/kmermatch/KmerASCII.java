/*
 * Copyright 2017 rad.
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
package kmermatch;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/**
 *
 * @author rad
 */
public class KmerASCII extends Kmer implements Comparable<KmerASCII>{
    private final byte[] ascii;

    public KmerASCII(String canonical) {
        ascii = canonical.getBytes(StandardCharsets.US_ASCII);
    }
    
    public byte[] getAscii() {
        return ascii;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(ascii);
    }

    @Override
    public boolean equals(Object obj) {
        final KmerASCII other = (KmerASCII) obj;
        return Arrays.equals(this.ascii, other.ascii);
    }
    
    @Override
    public int compareTo(KmerASCII another) {
        byte[] asciiOther = another.getAscii();
        int len = ascii.length < asciiOther.length ? ascii.length :asciiOther.length;
        for (int i = 0; i < len; i++) {
            if(ascii[i] < asciiOther[i]) {
                return -1;
            } else if(ascii[i] > asciiOther[i]) {
                return 1;
            }           
        }              
        return ascii.length - asciiOther.length;
    }
    
    
    
}
