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

import java.util.Collections;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class ParserTest {

    public ParserTest(String mainClassName, String moduleName, String[] args) {
        OptSet optSet = populateOptSet();
        UsageAndHelp usageAndHelp = optSet.getUsageAndHelp(mainClassName, moduleName, 120);
        System.err.println(usageAndHelp.getUsage());
        System.err.println(usageAndHelp.getHelp());
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, moduleName, 120);
        System.err.println("\n===================\noptions used and parsed\n===================");
        for (Opt o : optSet.getOptsList()) {
            if (o.isUsed()) {
                System.err.print(o.getOptLabelString() + "\t");
                for (Object x : o.getValues()) {
                    System.err.print(x + ", ");
                }
                System.err.println("");
            }
        }
        System.err.println("\n===================\nPositional options used and parsed\n===================");
        for (PositionalOpt o : optSet.getPositionalOptsList()) {
            if (o.isUsed()) {
                System.err.print(o.getUsageString() + "\t");
                if (o.getMaxOpts() == 1) {
                    System.err.println(o.getValue());
                } else {
                    for (Object x : o.getValues()) {
                        System.err.print(x + ", ");
                    }
                    System.err.println("");
                }
            }
        }

    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        optSet.setListingGroupLabel("[Input options]");
        optSet.addOpt(new Opt('a', "aaa", "Help for option a aaa"));
        optSet.addOpt(new Opt('b', null, "Help for option b"));
        optSet.addOpt(new Opt(null, "ccc", "Help for option ccc is very long ............ . ................. ....................... ................ .................... ........"));
        optSet.addOpt(new Opt('d', "d-very-very-very-long-key", "Help for option d"));
        optSet.addOpt(new Opt<>('e', "eeeeee", "Help for option e"));
        int currentGroup = optSet.incrementLisitngGroup();
        optSet.setListingGroupLabel(currentGroup, "[Output options]");
        optSet.addOpt(new Opt<>('x', "xxx", "Help for option x", 0.5, 0.0, 1.0, 1, 1));
        optSet.addOpt(new Opt<>('f', "fff", "Help for option f", 5, 0, 10, 1, 1));
        optSet.incrementLisitngGroup();
        Opt g = new Opt<>('G', "GGGGG", "Help for option G", 5, 0, 10, 3, 3);
        optSet.addOpt(g);
        boolean positionalArgumentRequired = true;
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME1", "name of input file ", 1, positionalArgumentRequired));
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME2", "name of input file ", 2));
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 3, 10));
        return optSet;
    }

}
