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

package shared;

/**
 *
 * @author rad
 */
public class Info {
        public static String getIupacCodesTable() {
        return ""+
"|---------+-------------------------------+----------------|\n" +
"|  Symbol | Description                   | Bases          |\n" +
"|---------+-------------------------------+---+---+---+----|\n" +
"|  A      | Adenine                       | A |   |   |    |\n" +
"|  C      | Cytosine                      |   | C |   |    |\n" +
"|  G      | Guanine                       |   |   | G |    |\n" +
"|  T      | Thymine                       |   |   |   | T  |\n" +
"|  U      | Uracil                        |   |   |   | U  |\n" +
"|  W      | Weak                          | A |   |   | T  |\n" +
"|  S      | Strong                        |   | C | G |    |\n" +
"|  M      | aMino                         | A | C |   |    |\n" +
"|  K      | Keto                          |   |   | G | T  |\n" +
"|  R      | puRine                        | A |   | G |    |\n" +
"|  Y      | pYrimidine                    |   | C |   | T  |\n" +
"|  B      | not A (B comes after A)       |   | C | G | T  |\n" +
"|  D      | not C (D comes after C)       | A |   | G | T  |\n" +
"|  H      | not G (H comes after G)       | A | C |   | T  |\n" +
"|  V      | not T (V comes after T and U) | A | C | G |    |\n" +
"|  N      | any Nucleotide                | A | C | G | T  |\n" +
"|---------+-------------------------------+---+---+---+----|\n" +
"Source: https://en.wikipedia.org/wiki/Nucleic_acid_notation\n";
        }
}
