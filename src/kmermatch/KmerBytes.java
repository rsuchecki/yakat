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
import kextender.CoreCoder;

/**
 *
 * @author rad
 */
public class KmerBytes extends Kmer implements Comparable<KmerBytes> {

    private final byte[] bytes;

    public KmerBytes(String canonical, boolean storeASCIII) {
        if (storeASCIII) {
            bytes = canonical.getBytes(StandardCharsets.US_ASCII);
        } else {
            bytes = CoreCoder.encodeCoreByteArray(canonical);
        }
//        int[] encodeCoreIntArray = CoreCoder.encodeCoreIntArray(canonical);
//        String decodeCore = CoreCoder.decodeCore(canonical.length(), encodeCoreIntArray);
//        System.err.println("IN : "+canonical);
//        System.err.println("OUT: "+decodeCore);
    }

    public byte[] getBytes() {
        return bytes;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(bytes);
    }

    @Override
    public boolean equals(Object obj) {
        final KmerBytes other = (KmerBytes) obj;
        return Arrays.equals(this.bytes, other.bytes);
    }

    @Override
    public int compareTo(KmerBytes another) {
        byte[] asciiOther = another.getBytes();
        return CoreCoder.compareCores(bytes, asciiOther);
//        int len = bytes.length < asciiOther.length ? bytes.length : asciiOther.length;
//        for (int i = 0; i < len; i++) {
//            if (bytes[i] < asciiOther[i]) {
//                return -1;
//            } else if (bytes[i] > asciiOther[i]) {
//                return 1;
//            }
//        }
//        return bytes.length - asciiOther.length; // N/A
    }

}
