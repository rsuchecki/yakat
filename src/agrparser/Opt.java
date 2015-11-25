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
package agrparser;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import shared.CommonMaths;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 * @param <T>
 */
public class Opt<T extends Comparable<T>> implements Comparable<Opt> {

    private T defaultValue;
    private T minValue;
    private T maxValue;
    private final Character shortKey;
    private final String longKey;
    private final String helpString;
    private Integer minValueArgs = 0;
    private int maxValueArgs = 0;
    private Integer listingGroup = 1;
    private Integer listingGroupPosition = 1;
    private boolean required;
    private boolean optFlag; //used to indicate that arg is used
    private HashMap<Integer, String> footnotesMap;

    //parsed values if opt is not a simple flag
    private final ArrayList<T> values = new ArrayList<>();

    /**
     * Create a "flag" opt
     *
     * @param shortKey [REQUIRED]
     * @param longKey [REQUIRED]
     * @param helpString [REQUIRED]
     */
    public Opt(Character shortKey, String longKey, String helpString) {
        this.shortKey = shortKey;
        this.longKey = longKey;
        this.helpString = helpString;
        
    }

    /**
     * Create an opt that takes maxRags args
     *
     * @param shortKey
     * @param longKey
     * @param helpString
     * @param numArgs
     */
    public Opt(Character shortKey, String longKey, String helpString, int numArgs) {
        this.shortKey = shortKey;
        this.longKey = longKey;
        this.helpString = helpString;
        this.minValueArgs = numArgs;
        this.maxValueArgs = numArgs;
    }

    /**
     * Constructor gives access to most fields, all but one of (shortKey,
     * longKey, helpString) may be null
     *
     * @param shortKey [REQUIRED]
     * @param longKey [REQUIRED]
     * @param helpString [REQUIRED]
     * @param defaultValue
     * @param minValue
     * @param maxValue
     * @param maxArgs [careful if >1 and positional args used]
     * @param listingGroup
     * @param listingGroupPosition
     */
    public Opt(Character shortKey, String longKey, String helpString, T defaultValue, T minValue, T maxValue,
            Integer maxArgs, Integer minArgs, Integer listingGroup, Integer listingGroupPosition) {
        this.shortKey = shortKey;
        this.longKey = longKey;
        this.helpString = helpString;
        this.defaultValue = defaultValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
        if (minArgs != null) {
            this.minValueArgs = minArgs;
        }
        if (maxArgs != null) {
            this.maxValueArgs = maxArgs;
        }
        if (listingGroup != null) {
            this.listingGroup = listingGroup;
        }
        if (listingGroupPosition != null) {
            this.listingGroupPosition = listingGroupPosition;
        }
    }

    /**
     * Constructor gives access to 3 required fields as well as the default and
     * min and max accepted values
     *
     * @param shortKey [REQUIRED]
     * @param longKey [REQUIRED]
     * @param helpString [REQUIRED]
     * @param defaultValue
     * @param minValue
     * @param maxValue
     */
    public Opt(Character shortKey, String longKey, String helpString, T defaultValue, T minValue, T maxValue) {
        this.shortKey = shortKey;
        this.longKey = longKey;
        this.helpString = helpString;
        this.defaultValue = defaultValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
        this.minValueArgs = 1;
        this.maxValueArgs = 1;
    }

    /**
     * Constructor gives access to 3 required fields as well as the default and
     * min and max accepted values plus the min max args
     *
     * @param shortKey [REQUIRED]
     * @param longKey [REQUIRED]
     * @param helpString [REQUIRED]
     * @param defaultValue
     * @param minValue
     * @param maxValue
     * @param minArgs
     * @param maxArgs [default=0]
     */
    public Opt(Character shortKey, String longKey, String helpString, T defaultValue, T minValue, T maxValue,
            Integer minArgs, Integer maxArgs) {
        this.shortKey = shortKey;
        this.longKey = longKey;
        this.helpString = helpString;
        this.defaultValue = defaultValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
        if (minArgs != null) {
            this.minValueArgs = minArgs;
        }
        if (maxArgs != null) {
            this.maxValueArgs = maxArgs;
        }
    }

    public void addValue(String valueString) {

        if (maxValueArgs == 0) {
            Reporter.reportNoMem("[FATAL]", "Boolean flag must not carry args. Offending value" + valueString, getClass().getSimpleName());
            System.exit(1);
        } else {        
            T value = (T) valueString;
            try {
                if (hasMinValue() && compareValues(getMinValue(), valueString) < 0) {
                    Reporter.reportNoMem("[FATAL]", "Value for option '" + getOptLabelString() + "' must be set to at least " + getMinValue() + ", offending value: " + valueString, getClass().getSimpleName());
                    System.exit(1);
                }
                if (hasMaxValue() && compareValues(getMaxValue(), valueString) > 0) {
                    Reporter.reportNoMem("[FATAL]", "Value for option '" + getOptLabelString() + "' must be set to at most " + getMaxValue() + ", offending value: " + valueString, getClass().getSimpleName());
                    System.exit(1);
                }
                if (getValues().size() >= maxValueArgs) {
                    Reporter.reportNoMem("[FATAL]", "Option '" + getOptLabelString() + "' allows up to " + getMaxValueArgs() + " values, unable to add : " + valueString, getClass().getSimpleName());
                    System.exit(1);
                }
                if(hasDefaultValue()) {
                    value = convertToType(getDefaultValue(), valueString);
                } else if (hasMinValue()){
                    value = convertToType(getMinValue(), valueString);                    
                } else if (hasMaxValue()){
                    value = convertToType(getMaxValue(), valueString);                    
                }
            } catch (ClassCastException | NumberFormatException e) {
                Reporter.reportNoMem("[FATAL]", "Possible reason(1): argument value type mismatch for option '" + getOptLabelString() + "', offending value: " + valueString, getClass().getSimpleName());
                Reporter.reportNoMem("[FATAL]", "Possible reason(2): a wrong number of values passed with option '" + getOptLabelString() + "', expected number of values: min=" + getMinValueArgs() + ", max=" + getMaxValueArgs() + "", getClass().getSimpleName());
                System.err.println(e.getMessage());
                System.exit(1);
            }
            values.add(value);
        }
    }
    
    private T convertToType(T storedValue, String inputValueString) throws NumberFormatException {
        Integer i = 0;
        Long l = 0L;
        Double d = 0.0;
        T value = (T) inputValueString;
        if (storedValue.getClass().isInstance(l)) {
            Long valueLong = Long.parseLong(inputValueString);
            value = (T) valueLong;
        } else if (storedValue.getClass().isInstance(i)) {
            Integer valueInteger = Integer.parseInt(inputValueString);
            value = (T) valueInteger;
        } else if (storedValue.getClass().isInstance(d)) {
            Double valueDouble = Double.parseDouble(inputValueString);
            value = (T) valueDouble;
        }
        return value;
    }

    private Integer compareValues(T storedValue, String inputValueString) throws NumberFormatException {
        T value = convertToType(storedValue, inputValueString);       
        return value.compareTo(storedValue);
    }

    public String getFormattedValue(T value) {
        Integer i = 0;
        Long l = 0L;
        Double d = 0.0;
        if (value.getClass().isInstance(d)) {
            d = (Double) value;
            if (Math.abs(d) < 10E9 && Math.abs(d) > 10E-9) {
                NumberFormat instance = DecimalFormat.getInstance();
                instance.setMaximumFractionDigits(9);
                return instance.format(d);
            } else {
                return CommonMaths.getScientific(d);
            }
        } else if (value.getClass().isInstance(i)) {
            i = (Integer) value;
            if (Math.abs(i) < 10E9 && Math.abs(i) > 10E-9) {
                return NumberFormat.getInstance().format(i);
            } else {
                return CommonMaths.getScientific(i);
            }
        } else if (value.getClass().isInstance(l)) {
            l = (Long) value;
            if (Math.abs(l) < 10E9 && Math.abs(l) > 10E-9) {
                return NumberFormat.getInstance().format(l);
            } else {
                return CommonMaths.getScientific(l);
            }
        } else {
            return (String) value;
        }

    }

    public Character getShortKey() {
        return shortKey;
    }

    public String getOptLabelString() {
        StringBuilder sb = new StringBuilder();
        sb.append(getOptLabelShort());
        if(hasShortKey() && hasLongKey()) {
            sb.append(" / ");
        }
        sb.append(getOptLabelLong());
        return sb.toString();
    }
    
    public String getOptLabelStringQuoted() {
        StringBuilder sb = new StringBuilder("'");
        sb.append(getOptLabelShort());
        if(hasShortKey() && hasLongKey()) {
            sb.append("' / '");
        }
        sb.append(getOptLabelLong()).append("'");
        return sb.toString();
    }
    
    public String getOptLabelShort() {
        StringBuilder sb = new StringBuilder();
        if (hasShortKey()) {
            sb.append("-").append(getShortKey());
            if (getMaxValueArgs() == 1) {
                sb.append(" <arg>");
            } else if (getMaxValueArgs() > 1) {
                sb.append(" <args>");
            }
        }
        return sb.toString();
    }
    
    public String getOptLabelLong() {
        StringBuilder sb = new StringBuilder();
        if (hasLongKey()) {
            sb.append("--").append(getLongKey());
            if (getMaxValueArgs() == 1) {
                sb.append(" <arg>");
            } else if (getMaxValueArgs() > 1) {
                sb.append(" <args>");
            }
        }
        return sb.toString();
    }

    public int getLongKeyLength() {
        if (longKey != null) {
            return longKey.length();
        } else {
            return 0;
        }
    }

    public String getLongKey() {
        return longKey;
    }

    public String getHelpString() {
        return helpString;
    }

    public T getDefaultValue() {
        return defaultValue;
    }

    public T getMinValue() {
        return minValue;
    }

    public T getMaxValue() {
        return maxValue;
    }

    public Integer getListingGroup() {
        return listingGroup;
    }

    public Integer getListingGroupPosition() {
        return listingGroupPosition;
    }

    public boolean hasShortKey() {
        return shortKey != null;
    }

    public boolean hasLongKey() {
        return longKey != null;
    }

    public boolean hasMinValue() {
        return minValue != null;
    }

    public boolean hasMaxValue() {
        return maxValue != null;
    }

    public boolean hasDefaultValue() {
        return defaultValue != null;
    }

    public boolean hasListingGroup() {
        return listingGroup != null;

    }

    public boolean hasListingGroupPosition() {
        return listingGroupPosition != null;
    }

    public Integer getMaxValueArgs() {
        return maxValueArgs;
    }

    public boolean takesArgs() {
        return getMaxValueArgs() > 0;
    }

    public Integer getMinValueArgs() {
        return minValueArgs;
    }

    public ArrayList<T> getValues() {
        return values;
    }

    public T getValueIfSingle() {
        if (values.size() == 1) {
            return values.get(0);
        }
        return null;
    }

    public T getValueOrDefault() {
        if (isUsed()) {
            if (values.size() == 1) {
                return values.get(0);
            }
        }        
        return getDefaultValue();
    }

    public boolean canTakeMoreValues() {
        return getValues().size() < getMaxValueArgs();
    }

    public void setListingGroup(Integer listingGroup) {
        this.listingGroup = listingGroup;
    }

    public void setListingGroupPosition(Integer listingGroupPosition) {
        this.listingGroupPosition = listingGroupPosition;
    }

    public boolean getOptFlag() {
        return optFlag;
    }

    public void setOptFlag(boolean optFlag) {
        this.optFlag = optFlag;
    }

    public boolean isUsed() {
        return optFlag || !getValues().isEmpty();
    }

    public int getNumberOfValues() {
        if (values == null) {
            return 0;
        } else {
            return values.size();
        }
    }

    @Override
    /**
     * Used for ordering opts in help output
     */
    public int compareTo(Opt another) {
        int interGroupOrder = getListingGroup() - another.getListingGroup();
        if (interGroupOrder == 0) {
            int intraGroupOrder = getListingGroupPosition() - another.getListingGroupPosition();
            if (intraGroupOrder == 0) {
                if (hasShortKey() && another.hasShortKey()) {
                    return getShortKey().compareTo(another.getShortKey());
                }
                if (hasLongKey() && another.hasLongKey()) {
                    return getLongKey().compareTo(another.getLongKey());
                }
            }
        }
        return interGroupOrder;
    }

    public boolean isVarArgsOption() {
        return getMinValueArgs().compareTo(getMaxValueArgs()) != 0;
    }

    public Opt setDefaultValue(T defaultValue) {
        this.defaultValue = defaultValue;
        return this;
    }

    public boolean isRequired() {
        return required;
    }

    public Opt setRequired(boolean required) {
        this.required = required;
        return this;
    }


    public Opt  setMinValue(T minValue) {
        this.minValue = minValue;
        return this;
    }

    public Opt  setMaxValue(T maxValue) {
        this.maxValue = maxValue;
        return this;
    }

    public Opt  setMinValueArgs(Integer minValueArgs) {
        this.minValueArgs = minValueArgs;
        return this;
    }

    public Opt  setMaxValueArgs(int maxValueArgs) {
        this.maxValueArgs = maxValueArgs;
        return this;
    }
    
    public Opt setNumArgs(int numArgs) {
        this.minValueArgs = numArgs;
        this.maxValueArgs = numArgs;
        return this;
    }


    public Opt addFootnote(int id, String footnoteText) {
        if(footnotesMap == null) {
            footnotesMap = new HashMap<>(2);
        }
        if(footnotesMap.containsKey(id)) {
            Reporter.reportNoMem("[FATAL]", "Footnote ["+id+"] for option '" + getOptLabelString() + " already set to "+footnotesMap.get(id), getClass().getSimpleName());
        } else {
            footnotesMap.put(id, footnoteText);
        }
        return this;
    }
    
    public ArrayList<Integer> getFootnoteKeys() {
        ArrayList<Integer> keys = new ArrayList<>(footnotesMap.size());
        for (Integer k: footnotesMap.keySet()) {
            keys.add(k);
        }
        Collections.sort(keys);
        return keys;
    }

    public HashMap<Integer, String> getFootnotesMap() {
        return footnotesMap;
    }
    
    public boolean hasFootnotes() {
        return footnotesMap != null;
    }
    
    public String getFootnote(int key) {
        return footnotesMap.get(key);
    }
   
    
}
