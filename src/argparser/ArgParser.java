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
import java.util.Iterator;
import java.util.ListIterator;
import shared.Reporter;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class ArgParser {

    /**
     * Supposedly POSIX-compliant CLI args processor
     *
     * @param args
     * @param optSet
     * @param firstArgModuleName
     * @param callerName
     * @param helpWidth
     */
    public void processArgs(String[] args, OptSet optSet, boolean firstArgModuleName, String callerName, int helpWidth) {
        boolean fatal = false;
        if (optSet.hasPositionalArgs() && optSet.hasVarArgOption()) {
            Reporter.reportNoMem("[FATAL]", "Mixing (i) options which take a variable number of values with (ii) positional arguments is asking for trouble", getClass().getSimpleName());
            fatal = true;
        }
        ArrayList<String> allArgs = new ArrayList<>();
        int argStart = 0;
        if (firstArgModuleName) {
            argStart++;
        }
        //SPLIT-UP POSIX STYLE (-abcd) IF ANY        
        for (int a = argStart; a < args.length; a++) {
            String arg = args[a];
            if (arg.startsWith("--")) {
                allArgs.add(arg);
            } else if (arg.startsWith("-")) {
                if (arg.length() > 2) {
                    char[] chars = arg.toCharArray();
                    StringBuilder sb = new StringBuilder();
                    for (int i = 1; i < chars.length; i++) {
                        char c = chars[i];
//                        System.out.println(c + " ASCII value " + (int) c);
                        if ((int) c >= 48 && (int) c <= 57) {
                            sb.append(c);
                        } else {
                            allArgs.add("-" + c);
                        }
                    }
                    String toString = sb.toString();
                    if (!toString.isEmpty()) {
                        allArgs.add(toString);
                    }

                } else {
                    allArgs.add(arg);
                }
            } else {
                allArgs.add(arg);
            }
        }

        //Now process args
        ListIterator<String> it = allArgs.listIterator();
        ArrayList<PositionalOpt> positionalOpts = optSet.getPositionalOptsList();
        ListIterator<PositionalOpt> positionalIterator = positionalOpts.listIterator();
        boolean parsingPositional = false;
        Opt opt = null;
        while (it.hasNext()) {
            String a = it.next();            
            if (a.equals("-h") || a.equals("--help")) {
                UsageAndHelp usageAndHelp = optSet.getUsageAndHelp(callerName, args[0], helpWidth);
                System.err.println();
                System.err.println(usageAndHelp.getUsage());
                System.err.println(usageAndHelp.getHelp());
                System.exit(0);
            }
            if (a.isEmpty()) {
                //skip empty fields
            } else if (a.startsWith("-") && !a.equals("-")) { //to allow dash stand-in for stdin: !a.equals("-")
                if (parsingPositional) {
                    Reporter.reportNoMem("[FATAL]", "Options mixed with positional arguments or incorrect number of values passed with option ", getClass().getSimpleName());
                    fatal = true;
                }
                Opt currentOpt = optSet.getOpt(a);
                if (currentOpt == null) {
                    System.err.println("Unrecognized argument " + a + ", try -h or --help. Terminating...");
                    fatal = true;
                }
                currentOpt.incrementOptInstance(); //accommodates for an opt that can be called multiple times 
                if (currentOpt.canTakeMoreValues()) {
                    opt = currentOpt;
                } else {
                    currentOpt.setOptFlag(true);
                }
            } else if (opt == null) {
                parsingPositional = true;
                if (positionalIterator.hasNext()) {
                    PositionalOpt positionalOpt = positionalIterator.next();
                    positionalOpt.addValue(a);
                    for (int i = 1; i < positionalOpt.getMaxOpts(); i++) {
                        if (it.hasNext()) {
                            String value = it.next();
                            positionalOpt.addValue(value);
                        } else {
                            break;
                        }
                    }
                } else {
                    Reporter.reportNoMem("[FATAL]", "Unexpected positional argument '" + a + "', try -h or --help", getClass().getSimpleName());
                    fatal = true;
                }
            } else {
                opt.addValue(a);
                if (!opt.canTakeMoreValues()) {
                    opt = null;
                }
            }
        }
        for (Opt o : optSet.getOptsList()) {
            if (o.isUsed() && o.getMinValueArgs() > o.getNumberOfValues()) {
                Reporter.reportNoMem("[FATAL]", "Insufficient (" + o.getNumberOfValues() + ") values passed with option '" + o.getOptLabelString() + "', at least " + o.getMinValueArgs() + " expected", getClass().getSimpleName());
                fatal = true;
            }
            if(o.isRequired() && !o.isUsed()) {
                Reporter.reportNoMem("[FATAL]", "Required option not used " + o.getOptLabelStringQuoted()+ " expected", getClass().getSimpleName());
                fatal = true;                
            }
        }
        if(fatal) {
            System.exit(1);
        }
    }
}
