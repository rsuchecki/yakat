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
package kextender;

/**
 *
 * @author rad
 */
public class SeedExtensionsPair {

    private String extensionLeft;
    private String extensionRight;

    public String getExtensionLeft() {
        return extensionLeft == null ? "" : extensionLeft;
    }

    public boolean setExtensionLeft(String extensionLeft) {
        if (this.extensionLeft == null || this.extensionLeft.length() < extensionLeft.length()) {
            this.extensionLeft = extensionLeft;
            return true;
        }
        return false;
    }

    public String getExtensionRight() {
        return extensionRight == null ? "" : extensionRight;
    }

    public boolean setExtensionRight(String extensionRight) {
        if (this.extensionRight == null || this.extensionRight.length() < extensionRight.length()) {
            this.extensionRight = extensionRight;
            return true;
        }
        return false;
    }

}
