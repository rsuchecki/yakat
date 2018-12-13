/*
 * Copyright 2017 Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au.
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
package allorfs;

import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.regex.Pattern;
import shared.Orf;
import shared.Reporter;
import shared.Sequence;
import shared.SequenceOps;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class OrfPredictorConsumer implements Runnable {

    private final BlockingQueue<ArrayList<Sequence>> inputQueue;
    private final String TOOL_NAME;
    private final int minLength;
    private final PrintStream bufferedOut;
    private final boolean translate;
    private final boolean requireStop;
    private final boolean strandMatcherPlus;
    private final boolean strandMatcherMinus;

    public OrfPredictorConsumer(BlockingQueue<ArrayList<Sequence>> inputQueue, String TOOL_NAME, int minLength, PrintStream bufferedOut,
            boolean translate, boolean requireStop, String strand) {
        this.inputQueue = inputQueue;
        this.TOOL_NAME = TOOL_NAME;
        this.minLength = minLength;
        this.bufferedOut = bufferedOut;
        this.translate = translate;
        this.requireStop = requireStop;
        strandMatcherPlus = Pattern.compile("both|plus|f(orward)?|fwd|\\+", Pattern.CASE_INSENSITIVE).matcher(strand).matches();
        strandMatcherMinus = Pattern.compile("both|minus|rev|r(everse)?|\\-", Pattern.CASE_INSENSITIVE).matcher(strand).matches();
        if(!strandMatcherPlus && !strandMatcherMinus) {
            Reporter.report("[FATAL]", "Ambiguous specpeciofication of -s / --strand by the user: " + strand, TOOL_NAME);
            System.exit(1);
        }        
    }

    @Override
    public void run() {
        ArrayList<Sequence> list;
        Pattern startStopCodonsForward = Pattern.compile("ATG|(T(AA|AG|GA))", Pattern.CASE_INSENSITIVE);
        Pattern startStopCodonsReverse = Pattern.compile("CAT|((TT|TC|CT)A)", Pattern.CASE_INSENSITIVE);
        Pattern stopCodonsForward = Pattern.compile("T(AA|AG|GA)", Pattern.CASE_INSENSITIVE);
        Pattern stopCodonsReverse = Pattern.compile("(TT|TC|CT)A", Pattern.CASE_INSENSITIVE);
        

                
        try {
            while (!(list = inputQueue.take()).isEmpty()) {
                for (Sequence sequence : list) {
                    Reporter.report("[INFO]", "Identifying ORFs on " + sequence.getId(), TOOL_NAME);
                    ArrayList<Orf> orfs = sequence.getOrfs(minLength, startStopCodonsForward, startStopCodonsReverse, stopCodonsForward, stopCodonsReverse);//, sequence, bufferedOut);
                    int count = 0;
                    for (Orf orf : orfs) {
//                System.out.printf("%8s%12d%12d%3d%12d\n",orf.getParenId(),orf.getFrom(),orf.getTo(),orf.getFrame(),orf.getLength());
                        if (!requireStop || orf.hasStopCodon()) {
                            if ((strandMatcherPlus && orf.getFrame() > 0) || (strandMatcherMinus && orf.getFrame() < 0)) {
                                CharSequence seq = sequence.getSequenceString().subSequence(orf.getFrom() - 1, orf.getTo());
                                StringBuilder sb = new StringBuilder(orf.getFastaHeader());
                                sb.append(System.lineSeparator());
                                if (orf.getFrame() < 0) {
                                    seq = SequenceOps.getReverseComplement(seq);
                                }
                                if (translate) {
                                    sb.append(SequenceOps.translate(seq));
                                } else {
                                    sb.append(seq);
                                }
                                int len = orf.getTo() - orf.getFrom() + 1;
                                int seqLength = seq.length();
                                if (seqLength != len) {
                                    Reporter.report("[BUG]", "ORF length mismatch", TOOL_NAME);
                                    System.err.println(sb);
                                    System.exit(1);
                                }
                                bufferedOut.println(sb);
                                count++;
                            }
                        }
                    }
                    Reporter.report("[INFO]", NumberFormat.getInstance().format(orfs.size()) + " ORFs found on " + sequence.getId(), TOOL_NAME);
                    Reporter.report("[INFO]", "Outputtig " + NumberFormat.getInstance().format(count) + " out of " + NumberFormat.getInstance().format(orfs.size()) + " ORFs identified on " + sequence.getId(), TOOL_NAME);
                }
            }
            inputQueue.put(new ArrayList<>()); //inform other threads
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
