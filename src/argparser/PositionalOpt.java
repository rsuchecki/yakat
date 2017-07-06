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
package argparser;

import java.util.ArrayList;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PositionalOpt implements Comparable<PositionalOpt> {

    private String usageString;
    private String helpString;
    private int position;
    private boolean required;
    private String value;
    private int maxOpts = 1;
    private ArrayList<String> values;

    public PositionalOpt(String usageString, String helpString, int position) {
        this.usageString = usageString;
        this.helpString = helpString;
        this.position = position;
    }

    /**
     * Only suitable for last positional opt if maxOpts > 1
     *
     * @param usageString
     * @param helpString
     * @param position
     * @param maxOpts
     */
    public PositionalOpt(String usageString, String helpString, int position, int maxOpts) {
        this.usageString = usageString;
        this.helpString = helpString;
        this.position = position;
        this.maxOpts = maxOpts;
    }

    public PositionalOpt(String usageString, String helpString, int position, boolean required) {
        this.usageString = usageString;
        this.helpString = helpString;
        this.position = position;
        this.required = required;
    }

    public PositionalOpt setRequired() {
        this.required = true;
        return this;
    }
       
    public boolean isRequired() {
        return required;
    }

    public String getUsageString() {
        return usageString;
    }

    public String getHelpString() {
        return helpString;
    }

    public int getPosition() {
        return position;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

    public void addValue(String value) {
        if (maxOpts > 1) {
            if (values == null) {
                values = new ArrayList<>(maxOpts);
            }
            values.add(value);
        } else {
            setValue(value);
        }
    }

    public int getMaxOpts() {
        return maxOpts;
    }

    public ArrayList<String> getValues() {
        return values;
    }
    
    public boolean isUsed() {
        return getValue() != null || (getValues() != null && !getValues().isEmpty());
    }

    @Override
    public int compareTo(PositionalOpt o) {
        return getPosition() - o.getPosition();
    }

}
