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

import argparser.OptSet;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 *
 * @author rad
 */
public class StdRedirect {

    public StdRedirect(OptSet optSet, String TOOL_NAME) {
        stdRedirect(optSet, TOOL_NAME);
    }
    
    private void stdRedirect(OptSet optSet, String TOOL_NAME) {
        String outRedirect;
        String errRedirect;
        if ((outRedirect = (String) optSet.getOpt("o").getValueOrDefault()) != null) {
            try {
                File file = new File(outRedirect);
                PrintStream printStream;
                printStream = new PrintStream(new FileOutputStream(file));
                System.setOut(printStream);
            } catch (FileNotFoundException ex) {
                Reporter.report("[ERROR]", "Failed redirecting stdout to " + outRedirect, TOOL_NAME);
            }
        }
        if ((errRedirect = (String) optSet.getOpt("e").getValueOrDefault()) != null) {
            try {
                File file = new File(errRedirect);
                PrintStream printStream;
                printStream = new PrintStream(new FileOutputStream(file));
                System.setErr(printStream);
            } catch (FileNotFoundException ex) {
                Reporter.report("[ERROR]", "Failed redirecting stderr to " + errRedirect, TOOL_NAME);
            }
        }
    }

}
