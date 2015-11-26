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

/**
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class Message {
    public enum Level {
        INFO, WARNING, ERROR, FATAL
    }
    private final Level level;
    private final String body;
    private final String caller;

    public Message(Level level, String body, String caller) {
        this.level = level;
        this.body = body;
        this.caller = caller;
    }

    public String getLevel() {
        return "["+level.toString()+"]";
    }

    public String getBody() {
        return body;
    }

    public String getCaller() {
        return caller;
    }

    
    
    
            
}
