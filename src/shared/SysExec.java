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
package shared;

import java.io.*;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class SysExec {


    public static String[] execute(String args) {
        String[] stdOutErr = new String[2];
        try {
            
            Process process = Runtime.getRuntime().exec(args.split(" "));
            BufferedReader stdOut = new BufferedReader(new InputStreamReader(process.getInputStream()));
            BufferedReader stdErr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
            process.waitFor();
            String line1 = null;
            String line2 = null;
            StringBuilder stdoutStringBuilder = new StringBuilder();
            while ((line1 = stdOut.readLine()) != null) {
                stdoutStringBuilder.append(line1).append(System.lineSeparator());
                System.out.println(line1);
            }
            stdOut.close();
            StringBuilder stderrStringBuilder = new StringBuilder();
            while ((line2 = stdErr.readLine()) != null) {
                stderrStringBuilder.append(line2).append(System.lineSeparator());
                System.err.println(line2);
            }            
            stdErr.close();
            stdOutErr[0] = stdoutStringBuilder.toString();
            stdOutErr[1] = stderrStringBuilder.toString();            
            System.err.println(process.exitValue());
            process.destroy();
        } catch (IOException | InterruptedException e) {
            System.err.println(e.getMessage());
        }
        System.exit(1);
        return stdOutErr;
    }

}
