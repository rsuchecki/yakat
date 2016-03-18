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

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import shared.CommonMaths;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class OptSet {

    HashMap<Character, Opt> shortToOptMap;
    HashMap<String, Opt> longToOptMap;
    HashMap<Integer, ArrayList<Opt>> groupToOptsMap;
    ArrayList<Opt> optsList;

    HashMap<Integer, PositionalOpt> positionalOptsMap;
    ArrayList<PositionalOpt> positionalOptsList;

    Integer currentListingGroup = 0;
    Integer currentListingGroupPosition = 0;

    HashMap<Integer, String> listingGroupLabels;     
    HashMap<Integer, ArrayList<Integer>> mutuallyExclusiveGroups;
    ArrayList<Integer> requiredGroups;
    
    String summary;
    

//    ArrayList<OptGroup> optGroups; //taking care of mutual exlusion of groups of opts
    
    
//    /**
//     * Add an "flag" opt
//     * @param shortKey
//     * @param longKey
//     * @param helpString
//     * @return 
//     */
//    public boolean addOpt(Character shortKey, String longKey, String helpString) {
//        Opt o = new Opt(shortKey, longKey, helpString);
//        return addOpt(o);
//    }
//    
    public OptSet() {
        shortToOptMap = new HashMap<>();
        longToOptMap = new HashMap<>();
        optsList = new ArrayList<>();
        groupToOptsMap = new HashMap<>();
        positionalOptsMap = new HashMap<>();
        positionalOptsList = new ArrayList<>();
        requiredGroups = new ArrayList<>();
//        optGroups = new ArrayList<>();
    }
    
    public OptSet(String summary) {
        shortToOptMap = new HashMap<>();
        longToOptMap = new HashMap<>();
        optsList = new ArrayList<>();
        groupToOptsMap = new HashMap<>();
        positionalOptsMap = new HashMap<>();
        positionalOptsList = new ArrayList<>();
        requiredGroups = new ArrayList<>();
        this.summary = summary;
//        optGroups = new ArrayList<>();
    }

    public boolean hasVarArgOption() {
        for (Opt o : getOptsList()) {
            if (o.isVarArgsOption()) {
                return true;
            }
        }
        return false;
    }

    public boolean hasPositionalArgs() {
        return !positionalOptsList.isEmpty();
    }

    public void addPositionalOpt(PositionalOpt p) {
        if (positionalOptsMap.containsKey(p.getPosition())) {
            System.err.println("Failed addPositionalOpt() operation for " + p.getUsageString());
            System.err.println("Another positional argument already set for position" + p.getPosition());
            System.exit(1);
        } else {
            positionalOptsMap.put(p.getPosition(), p);
            positionalOptsList.add(p);
        }
    }

//    public void addOptGroup(OptGroup optGroup) {
//        optGroups.add(optGroup);
//        for (Opt opt : optGroup.getOpts()) {
//            addOpt(opt);
//        }
//    }
    
    public void addOpt(Opt opt) {
        addOpt(opt, currentListingGroup, currentListingGroupPosition++);
    }

    public void addOpt(Opt opt, Integer listingGroup, Integer listingGroupPosition) {
        Character shortKey = opt.getShortKey();
        String longKey = opt.getLongKey();
        boolean success = true;
        if (shortKey == null && longKey == null) {
//            return false;
            success = false;
        }

        if (shortKey != null) {
            if (shortToOptMap.containsKey(shortKey)) {
                success = false;
//            return false;
            } else {
                shortToOptMap.put(shortKey, opt);
            }
        }
        if (longKey != null) {
            if (longToOptMap.containsKey(longKey)) {
                success = false;
//            return false;
            } else {
                longToOptMap.put(longKey.toLowerCase(), opt);
            }
        }
        if (!success) {
            System.err.println("Failed addOpt() operation for " + opt.getOptLabelString() + " due to arg key clash");
            System.exit(1);
        }
        if (listingGroup != null) {
            opt.setListingGroup(listingGroup);
        }
        if (listingGroupPosition != null) {
            opt.setListingGroupPosition(listingGroupPosition);
        }
        optsList.add(opt);
//        return true;
    }

    public Opt getOptS(Character shortArg) {
        return shortToOptMap.get(shortArg);
    }

    public Opt getOptL(String longArg) {
        return longToOptMap.get(longArg.toLowerCase());
    }

    public ArrayList<Opt> getOptsList() {
        return optsList;
    }

    public ArrayList<PositionalOpt> getPositionalOptsList() {
        Collections.sort(positionalOptsList);
        for (int i = 0; i < positionalOptsList.size() - 1; i++) {
            PositionalOpt pos = positionalOptsList.get(i);
            if (pos.getMaxOpts() > 1) {
                System.err.println("Incoherent positional args list. Only the last positional option may take multiple inputs");
                System.exit(1);
            }
        }
        return positionalOptsList;
    }

    public Opt getOpt(String key) {
        Opt opt = null;
        while (key.startsWith("-")) {
            key = key.replaceFirst("-", "");
        }
        if (key.length() == 1) {
            opt = shortToOptMap.get(key.charAt(0));
        } else if (key.length() > 1) {
            opt = longToOptMap.get(key.toLowerCase());
        }
        if (opt == null) {
            Reporter.report("[FATAL]", "Unknown option requested: " + key, this.getClass().getSimpleName());
            System.exit(1);
        }
        return opt;
    }

    /**
     *
     * @param mainClassName
     * @param moduleName
     * @param printWidth
     * @return
     */
    public UsageAndHelp getUsageAndHelp(String mainClassName, String moduleName, int printWidth) {
        Collections.sort(getOptsList());
        Collections.sort(getPositionalOptsList());
        int maxLongArgLength = 0;
        for (Opt opt : getOptsList()) {
            if (opt.hasLongKey()) {
                maxLongArgLength = Math.max(maxLongArgLength, opt.getLongKey().length());
            }
        }
        maxLongArgLength++;
        int offset = 8 + maxLongArgLength + 8;

//        Map<String, String> env = System.getenv();
//        for (String envName : env.keySet()) {
//             System.out.format("%s=%s%n", envName, env.get(envName));
//        }
        String terminalColumns = System.getenv("COLUMNS");
        if (terminalColumns != null) {
            printWidth = Integer.parseInt(terminalColumns);
        }

        int helpLineWidth = printWidth - offset;

        //Generate usage string
        StringBuilder usage = new StringBuilder();
        if(summary != null) {
            usage.append(Reporter.wrapString(summary, printWidth)).append(System.lineSeparator());
        }
        usage.append(System.lineSeparator()).append("java -jar ").append(mainClassName).append(".jar ").append(moduleName);

        //Generate help page
        StringBuilder help = new StringBuilder();
        int currentListingGroup = -1;
        boolean hasRequiredOption = false;

        HashMap<Integer, String> footnotes = new HashMap<>();
        boolean wrappedLine = false;
        int optInUsage = 0;
        for (int i = 0; i < getOptsList().size(); i++) {
            Opt opt = getOptsList().get(i);
            //Split option groups
            if (opt.getListingGroup() > currentListingGroup) {
                currentListingGroup = opt.getListingGroup();
                String groupLabel = getListingGroupLabel(currentListingGroup);
                if (groupLabel.isEmpty()) {
                    if (i != 0) {
                        help.append(String.format("%" + (offset - 1) + "s", " ")).append(": ").append(System.lineSeparator());
                    }
                } else {
                    //if last line not wrapped add extra line for better readability
                    if (!wrappedLine) {
                        help.append(System.lineSeparator());
                    }
                    help.append(groupLabel); 
                    ArrayList<Integer> mutuallyExclusive = getMutuallyExclusive(currentListingGroup);
                    if(isGroupRequired(currentListingGroup)) {
                        help.append("[**] ");
                    }
//                    if(mutuallyExclusive !=null && !mutuallyExclusive.isEmpty()) {
//                        help.append(" Mutually exclusive with:");
//                        for (Integer g : mutuallyExclusive) {
//                            help.append(" ").append(getListingGroupLabel(g).replaceAll("\\[|\\]", ""));
//                        }
//                    }
                    help.append(System.lineSeparator());
                }
            }
            help.append(" ");
            if (opt.hasShortKey()) {
                help.append("-").append(opt.getShortKey());
            } else {
                help.append("    ");
            }
            if (opt.hasShortKey() && opt.hasLongKey()) {
                help.append(", ");
            }
            if (opt.hasLongKey()) {
                help.append("--").append(opt.getLongKey());
            } else {
//                help.append("    ");
            }
            if (opt.getMaxValueArgs() > 1 || opt.getMinValueArgs() > 1) {
                help.append(" <args>");
            } else if (opt.getMaxValueArgs() == 1 || opt.getMinValueArgs() == 1) {
                help.append(" <arg> ");
            } else {
                help.append("       ");
            }
            if (!opt.hasLongKey()) {
                help.append("    ");                
            }

            String gap = String.format("%" + (maxLongArgLength - opt.getLongKeyLength()) + "s", " ");
            help.append(gap).append(": ");
            StringBuilder helpLine = new StringBuilder();
            if (opt.isRequired()) {
                helpLine.append("[*] ");
                hasRequiredOption = true;
                if (opt.hasShortKey()) {
                    usage.append(" ").append(opt.getOptLabelShort());
                    optInUsage++;
                }
            }
            helpLine.append(opt.getHelpString());
            if (opt.hasMinValue()) {
                if (helpLine.length() > 0) {
                    helpLine.append("; ");
                }
                helpLine.append("min=").append(opt.getFormattedValue(opt.getMinValue()));
            }
            if (opt.hasMaxValue()) {
                if (helpLine.length() > 0) {
                    helpLine.append("; ");
                }
                helpLine.append("max=").append(opt.getFormattedValue(opt.getMaxValue()));
            }
            if (opt.hasDefaultValue()) {
                if (helpLine.length() > 0) {
                    helpLine.append("; ");
                }
                helpLine.append("default=").append(opt.getFormattedValue(opt.getDefaultValue()));
            } else if (opt.getMaxValueArgs() > 0 && !opt.isRequired()) {
                if (helpLine.length() > 0) {
                    helpLine.append("; ");
                }
                helpLine.append("no default value");
            }
            if (opt.hasFootnotes()) {
                helpLine.append("; ");
                ArrayList<Integer> keys = opt.getFootnoteKeys();
                for (Integer key : keys) {
                    helpLine.append("[").append(key).append("]");
                    footnotes.put(key, opt.getFootnote(key));
                }
            }
            help.append(Reporter.wrapString(helpLine.toString(), helpLineWidth, offset));
            //record if line was wrapped
            if (helpLine.length() > helpLineWidth - offset) {
                wrappedLine = true;
            } else {
                wrappedLine = false;
            }
            help.append(System.lineSeparator());
        }
        if (hasRequiredOption) {
            help.append(System.lineSeparator()).append("[*] Denotes a required setting/argument.");
        }
        //Print footnotes
        if (!footnotes.isEmpty()) {
            ArrayList<Integer> keys = new ArrayList<>(footnotes.size());
            for (Integer k : footnotes.keySet()) {
                keys.add(k);
            }
            Collections.sort(keys);
            for (Integer key : keys) {
                help.append(System.lineSeparator()).append("[").append(key).append("] ").append(footnotes.get(key));
            }
        }

        //finish building usage string
        if (getOptsList().size() > optInUsage) {
            usage.append(" [options]");
        }
        for (PositionalOpt positional : getPositionalOptsList()) {
            if (positional.isRequired()) {
                usage.append(" ").append(positional.getUsageString());
            } else {
                usage.append(" [").append(positional.getUsageString()).append("]");
            }
        }
        String usageString = usage.toString();
        String helpString = help.toString().replaceAll("\\. *,", ",").replaceAll(" *,", ",").replaceAll(" *;", ";");
        
        UsageAndHelp usageAndHelp = new UsageAndHelp(usageString, helpString);
        return usageAndHelp;
        //        String s = "Currently k-mer frequency is not taken into consideration, so use of a dedicated k-mer counting program, "
//                + "such as KMC or Jellyfish is recommended. It is best to exclude low frequency k-mers before passing "
//                + "the list of k-mers to KmerExtender. For smaller jobs FASTA or FASTQ input may suffice.";
//        System.out.println(Reporter.wrapString(s, 145));
    }

    /**
     * Increment and
     *
     * @return current listing group number
     */
    public int incrementLisitngGroup() {
        return ++currentListingGroup;
    }

    public void setListingGroupLabel(String label) {
        setListingGroupLabel(currentListingGroup, label);
    }

    /**
     *
     * @param group
     * @param label (e.g. "Input options")
     */
    public void setListingGroupLabel(Integer group, String label) {
        if (listingGroupLabels == null) {
            listingGroupLabels = new HashMap<>();
        }
        listingGroupLabels.put(group, label);
    }

    public String getListingGroupLabel(Integer group) {
        if (listingGroupLabels != null && listingGroupLabels.get(group) != null) {
            return listingGroupLabels.get(group);
        }
        return "";
    }
    
    
    public void setGroupRequired(int group) {
        requiredGroups.add(group);
    }
    
    public boolean isGroupRequired(int group) {
        return requiredGroups.contains(group);
    }
    
    public void setMutuallyExclusiveGroups(int group1, int group2) {
        if(mutuallyExclusiveGroups == null) {
            mutuallyExclusiveGroups = new HashMap<>();                        
        }        
        setMustuallyExclusive(group1, group2);
        setMustuallyExclusive(group2, group1);
    }
    
    private void setMustuallyExclusive(int group1, int group2) {
        ArrayList<Integer> incompatibleWithGroup1 = mutuallyExclusiveGroups.get(group1);
        if(incompatibleWithGroup1 == null) {
            incompatibleWithGroup1 = new ArrayList<>();
        }
        incompatibleWithGroup1.add(group2);
        mutuallyExclusiveGroups.put(group1, incompatibleWithGroup1);
    }
    
    private  ArrayList<Integer> getMutuallyExclusive(int group) {
        if(mutuallyExclusiveGroups == null) {
            return null;
        }
        return mutuallyExclusiveGroups.get(group);
    }

}
