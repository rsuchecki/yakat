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
package gbssplit;

import shared.SequenceOps;

/**
 *
 * @author rad
 */
public class Sample {
    private final String barcode;
    private final String barcodeRC;
    private final String barcodeAndRestrictionRemnant1;
    private final String restrictionRemnant1;
    private final String restrictionRemnant1RC;
    private final String sampleId;

    public Sample(String barcode, String restrictionRemant1, String sampleId) {
        this.barcode = barcode;
        this.barcodeRC = SequenceOps.getReverseComplementString(barcode);
        this.barcodeAndRestrictionRemnant1 = barcode + restrictionRemant1;
//        this.batcodeAndRestrictionRemnant1 = barcode + "TGCAG";
        this.sampleId = sampleId;
        this.restrictionRemnant1 = restrictionRemant1;
        this.restrictionRemnant1RC = SequenceOps.getReverseComplementString(restrictionRemant1);
        
    }

    public String getBarcode() {
        return barcode;
    }

    public String getBarcodeRC() {
        return barcodeRC;
    }

    public String getRestrictionRemnant1() {
        return restrictionRemnant1;
    }

    public String getRestrictionRemnant1RC() {
        return restrictionRemnant1RC;
    }

    
    public String getBarcodeAndRestrictionRemnant1() {
        return barcodeAndRestrictionRemnant1;
    }

    public String getSampleId() {
        return sampleId;
    }
    
    
    
}
