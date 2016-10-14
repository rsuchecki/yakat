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
    private final String batcodeAndPstI;
    private final String sampleId;

    public Sample(String barcode, String sampleId) {
        this.barcode = barcode;
        this.barcodeRC = SequenceOps.getReverseComplementString(barcode);
        this.batcodeAndPstI = barcode + "TGCAG";
        this.sampleId = sampleId;
    }

    public String getBarcode() {
        return barcode;
    }

    public String getBarcodeRC() {
        return barcodeRC;
    }

    public String getBatcodeAndPstI() {
        return batcodeAndPstI;
    }

    public String getSampleId() {
        return sampleId;
    }
    
    
    
}
