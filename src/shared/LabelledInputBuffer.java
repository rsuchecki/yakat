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

import java.util.ArrayList;

/**
 *
 * @author rad
 */
public class LabelledInputBuffer {
    private final  String label;
    private final ArrayList data;

    public LabelledInputBuffer(String label, ArrayList data) {
        this.label = label;
        this.data = data;
    }

    public String getLabel() {
        return label;
    }

    public ArrayList getData() {
        return data;
    }
    
    
}
