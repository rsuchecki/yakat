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

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class CommonMaths implements Serializable {

    public static double round(double value, int decimalPoints) {
        double scaleFactor = Math.pow(10, decimalPoints);
        double scaledUp = value * scaleFactor;
        double rounded = Math.round(scaledUp);
        double scaledDown = rounded / scaleFactor;
//        System.out.println("\nAveraged = "+scaledDown);
        return scaledDown;
    }
    
    public static int round(double value) {
        return (int) round(value, 0);
    }
    

    public static String getRandomString() {
        StringBuilder sb = new StringBuilder("");
        Random generator = new Random();
        for (int i = 0; i < 4; i++) {
            sb.append((char) (48 + generator.nextInt(10)));
            sb.append((char) (65 + generator.nextInt(26)));
        }
        return sb.toString();

    }

//    1 Byte = 8 Bit
//1 Kilobyte = 1,024 Bytes
//1 Megabyte = 1,048,576 Bytes
//1 Gigabyte = 1,073,741,824 Bytes
//1 Terabyte = 1,099,511,627,776 Bytes
    public static String getBytesMultiple(long bytes) {
        if (bytes < 1024) {
            return bytes + " Bytes";
        } else if (bytes < 1048576) {
            return NumberFormat.getInstance().format(CommonMaths.round((double) bytes / 1024, 2)) + " kB";
        } else if (bytes < 1073741824) {
            return NumberFormat.getInstance().format(CommonMaths.round((double) bytes / 1048576, 2)) + " MB";
        } else {
            return NumberFormat.getInstance().format(CommonMaths.round((double) bytes / 1073741824, 2)) + " GB";
        }
    }

    public static String getScientific(Object number) {
        DecimalFormat df = new DecimalFormat("0.00E00");
        return df.format(number).toLowerCase();
    }
    
    public static double getMedian(ArrayList<Short> nonZeroValues) {
        if (nonZeroValues.isEmpty()) {
            return 0;
        } else {
            Collections.sort(nonZeroValues);
            int size = nonZeroValues.size();
            if (size == 1) {
                return (double) nonZeroValues.get(0);
            } else if (size % 2 == 0) {
                int middle = size / 2;
                int sum = nonZeroValues.get(middle - 1) + nonZeroValues.get(middle);
                return ((double) (sum)) / 2;
            } else {
                int middle = size / 2;
                return (double) nonZeroValues.get(middle);
            }
        }
    }
}
