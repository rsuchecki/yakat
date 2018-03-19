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
public class AuxiliaryInfoCoder {
    //TO ENCODE 
    //countLeftA
    //countLeftC
    //countLeftG
    //countLeftT
    //countRightA
    //countRightC
    //countRightG
    //countRightT
    
    
    
    
    public static long encode(Character c, boolean isLeft, int frequency) {
        
//    private char clipLeft = '#';   //2B                     ///*TODO*/: encode to 2-3b if sticking to int array
//    private char clipRight = '#';  //2B                     //TODO: encode to 2-3b
////    private byte storedCount;  //1B    //can switch to short since all those add up to 7 bytes, 1 wasted
//    private byte storedCountLeft;
//    private byte storedCountRigth;
//    private boolean invalid;  //1B
////    private boolean visited;//1B    
//    //then round to multi of 8
//
//    private byte visitedBy = Byte.MAX_VALUE;
return 0L;
    } 
}
