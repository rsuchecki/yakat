/*
 * Copyright 2016 Australian Centre For Plant Functional Genomics
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

package argparser;

import java.util.ArrayList;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class OptGroup {
    ArrayList<OptGroup> incompatibleOptGroups;
    ArrayList<Opt> opts;

    public OptGroup() {
        opts = new ArrayList<>();
    }   
    
    public OptGroup(ArrayList<Opt> opts) {
        this.opts = opts;
    }
    
    public void addOpt(Opt opt) {
        opts.add(opt);
    }

    public void addIncompatibleOptGroup(OptGroup incompatibleOptGroup) {
        if(incompatibleOptGroups == null) {
            incompatibleOptGroups = new ArrayList<>();
        }
        this.incompatibleOptGroups.add(incompatibleOptGroup);
    }

    public ArrayList<OptGroup> getIncompatibleOptGroups() {
        return incompatibleOptGroups;
    }

    public ArrayList<Opt> getOpts() {
        return opts;
    }
    
    
   
}
