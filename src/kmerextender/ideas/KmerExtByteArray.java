/*
 * Copyright 2015 Australian Centre For Plant Functional Genomics
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
package kmerextender.ideas;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class KmerExtByteArray {

    private final byte[] array; 

    public KmerExtByteArray(int length) {
        array = new byte[length];
    }

    public void add(byte encodedMer, int coreHashCode) {
        //CONVERT HASH OF THE PairMer CORE TO INDEX IN OUR ARRAY
        int index = coreHashCode % array.length;
        if (index < 0) {
            index = -index;
        }
        //STORE THE NEW VALUE IF ARR LOCATION EMPTY
        addSynchronized(encodedMer, index);
        

    }

    private synchronized void addSynchronized(byte encodedMer, int index) {
        if(array[index] == 0) {
            array[index] = encodedMer;
        } else {
            
        }
        
    }
    
    
    
   

}
