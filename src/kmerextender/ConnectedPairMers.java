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
package kmerextender;

import java.util.ArrayList;
import shared.SequenceOps;
import shared.Reporter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class for storing a set of connected PairMers with the connections between
 * them
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class ConnectedPairMers {
//    private ArrayList<PairMerNode> connectedPairMers = new ArrayList<>();

    private PairMerNode terminal1;
    private PairMerNode terminal2;
    private PairMerNode singletonNode;
    private final HashMap<PairMer, PairMerNode> pairMerNodes = new HashMap<>();
//    ArrayList<PairMer> connectedPairMersTemp = new ArrayList<>();
    private CharSequence extendedString;
    private boolean potentialDuplicate = false;
    

    /**
     * Given a PairMer recursively follows the implicit connections to generate
     * a slightly more explicit graph structure, which can than be quickly
     * traversed to generate the extended string using toString(int k)
     *
     * @param pairMer
     * @param k
     * @param pairMersMap
     * @param threadId
     * @return
     */
    public boolean connectPairMers(PairMer pairMer, int k, PairMersMap pairMersMap, byte threadId, String DEBUG_FILE) throws StackOverflowError {
//        boolean visited;
//        boolean visited = pairMer.checkAndSetVisited();
//        if (!pairMer.isVisited()) {
//            pairMer.setVisited();
        byte visitedBy = pairMer.checkAndSetVisitedBy(threadId);
        //????
        if (visitedBy == Byte.MAX_VALUE) { //NODE NOT VISITED BEFORE -> keep going
//            visited = false;
        } else if (visitedBy > threadId) { //NODE PREVIOUSLY VISITED BY A LOWER PRIORITY THREAD -> keep going, but tag result as a posiible duplicate
//            visited = false;
            potentialDuplicate = true;
        } else if (visitedBy == threadId) { //NODE ALREADY VISITED BY THE SAME THREAD -> skip
//            visited = true;
            return true;
        } else { //NODE PREVIOUSLY VISITED BY A HIGHER PRIORITY THREAD  -> terminate extension 
            return false;
        }

//        if(checkAndSetVisitedBy)
//        if (!visited || reset) {
//        pairMer.setVisited();
//            pairMer.setVisitedBy(threadId);

        String decodedCore = pairMer.decodeCore(k - 1);
        StringBuilder otherCoreOfKmer1 = new StringBuilder();
        StringBuilder otherCoreOfKmer2 = new StringBuilder();
        if (pairMer.hasLeftClip()) {
            otherCoreOfKmer1.append(pairMer.getClipLeft());
            otherCoreOfKmer1.append(decodedCore.subSequence(0, decodedCore.length() - 1));
        }
        if (pairMer.hasRightClip()) {
            otherCoreOfKmer2.append(decodedCore.subSequence(1, decodedCore.length()));
            otherCoreOfKmer2.append(pairMer.getClipRight());
        }
        PairMer otherPairMer1 = null;
        PairMer otherPairMer2 = null;

//        if(otherCoreOfKmer1.length() ==0 || otherCoreOfKmer2.length() ==0) {
//            int x =0;
//        }
        try {
            otherPairMer1 = pairMersMap.get(otherCoreOfKmer1, k);
            otherPairMer2 = pairMersMap.get(otherCoreOfKmer2, k);
        } catch (NonACGTException ex) {
            Reporter.report("[WARNING]", "Unexpected NonACGTException caught", getClass().getCanonicalName());
        }

        //SWITCH TO bit-level RC
//            PairMer otherPairMer1 = pairMer.getOtherPairmerCoreLeft(k);
//            PairMer otherPairMer2 = pairMer.getOtherPairmerCoreRight(k);
//            
//                otherCoreRight = pairMer.shiftBitsLeftAndAddRightClip(pairMer, k-1);
//            String[] otherCores = pairMer.getOtherCores(k);
//            PairMer tempLeft = pairMer.getOtherPairmerCoreLeft(k);
//            PairMer tempRight = pairMer.getOtherPairmerCoreRight(k);
//            if (otherPairMer1 != null && otherPairMer2 != null) {
//                if(!tempLeft.equals(otherPairMer1) || !tempRight.equals(otherPairMer2)) {
//                    System.err.println(otherCoreOfKmer1+" <- otherCore1");
//                    System.err.println(otherPairMer1.decodeCore(k - 1)+" <- via String (1)");                    
//                    System.err.println(tempLeft.decodeCore(k - 1)+" <- encoded");
//                    System.err.println("");
//                    System.err.println(otherCoreOfKmer2+" <- otherCore2");
//                    System.err.println(otherPairMer2.decodeCore(k - 1)+" <- via String (2)");                    
//                    System.err.println(tempRight.decodeCore(k - 1)+" <- encoded");
//                    System.err.println("");
//                    System.err.println("--------------------------------");
//                }
////            System.err.println(otherCoreOfKmer2);
//            }
//            if(pairMer.getPairMerString(k).equals("TACAAATACATATCCTTAACATACAAGATCAATGATAGAGAACGTG")) {
//            System.err.println(pairMer.getPairMerString(k) + "\tPairMer");
//            System.err.println(wrapper.getKmer1String() + "\tk1");
//            System.err.println(SequenceOps.getReverseComplementString(wrapper.getKmer1String()) + "\tk1 RC");
////            NavigableSet<PairMer> keySet = pairMersMap.getPairMersSkipListMap().keySet();
////            for (PairMer pairMer1 : keySet) {
////                System.err.println(pairMer1.getClipLeft()+"_"+pairMer1.getTmpCore()+"_"+pairMer1.getClipRight());
////            }
//            }
        boolean otherPairMer1isRC = false;
        if (otherPairMer1 != null) {

//                System.err.println(wrapper.getOtherCoreOfKmer1() + "\tk1otherCore");
//                System.err.println(SequenceOps.getReverseComplementString(wrapper.getOtherCoreOfKmer1()) + "\tk1otherCoreRC");
//                if (!otherCoreOfKmer1.equals(otherPairMer1.decodeCore(k - 1))) {
            if (!otherCoreOfKmer1.toString().equals(otherPairMer1.decodeCore(k - 1))) {
                otherPairMer1isRC = true;
//                    System.err.println("blah!!!!");
            }

//                //recycling variables visited, visitedBy
//                visitedBy = otherPairMer1.checkAndSetVisitedBy(threadId);
//                if (visitedBy == Byte.MAX_VALUE) { //NODE NOT VISITED BEFORE -> keep going
//                    visited = false;
//                } else if (visitedBy > threadId) { //NODE PREVIOUSLY VISITED BY A LOWER PRIORITY THREAD -> keep going, but tag result as a posiible duplicate
//                    visited = false;
//                    potentialDuplicate = true;
//                } else if (visitedBy == threadId) { //NODE ALREADY VISITED BY THE SAME THREAD -> skip
//                    visited = true;
//                } else { //NODE PREVIOUSLY VISITED BY A HIGHER PRIORITY THREAD  -> terminate extension 
//                    return false;
//                }
//                if (!visited) {
            if (!connectPairMers(otherPairMer1, k, pairMersMap, threadId, DEBUG_FILE)) {
                return false;
            }
//                    terminalPairMer1 = addPairMerToListOfConnected(otherPairMer1, k, pairMersMap, connectedPairMers);
//                }
        }
        boolean otherPairMer2isRC = false;
//            System.err.println(wrapper.getKmer2String() + "\tk2");
//            System.err.println(SequenceOps.getReverseComplementString(wrapper.getKmer2String()) + "\tk2 RC");
        if (otherPairMer2 != null) {
//                System.err.println(wrapper.getOtherCoreOfKmer2() + "\tk2otherCore");
//                System.err.println(SequenceOps.getReverseComplementString(wrapper.getOtherCoreOfKmer2()) + "\tk2otherCoreRC");
            if (!otherCoreOfKmer2.toString().equals(otherPairMer2.decodeCore(k - 1))) {
                otherPairMer2isRC = true;
            }
//                if (!otherPairMer2.isVisited()) {
            if (!connectPairMers(otherPairMer2, k, pairMersMap, threadId, DEBUG_FILE)) {
                return false;
            }
//                }
        }
        add(pairMer, otherPairMer1, otherPairMer1isRC, otherPairMer2, otherPairMer2isRC, k, DEBUG_FILE);
//        }
        return true;
    }

//    public boolean isBug() {
//        return bug;
//    }
    /**
     *
     * @param pairMer
     * @param previous
     * @param previousRc
     * @param next
     * @param nextRc
     */
    private void add(PairMer pairMer, PairMer previous, boolean previousRc, PairMer next, boolean nextRc, int k, String DEBUG_FILE) {
        PairMerNode pairMerNode = new PairMerNode(pairMer, previous, previousRc, next, nextRc);
//        connectedPairMers.add(pairMerNode);
//        if (previous == null ^ next == null) {
        if ((previous == null && next != null) || (previous != null && next == null)) {
            if (terminal1 == null) {
                terminal1 = pairMerNode;                
            } else if (terminal2 == null) {
                terminal2 = pairMerNode;
            } else {
                Reporter.report("[BUG?]", "Third terminal PairMerNode identified? Trying to store a terminal node again?", getClass().getSimpleName());
//                bug = true;
                if (DEBUG_FILE != null) {
                    ArrayList<String> toReport = new ArrayList<>();
                    try {

                        toReport.add("Third terminal PairMerNode identified? Trying to store a terminal node again?");

//                    System.err.println("Stored connected pairMers:");
//                    for (PairMer pairMerTmp : connectedPairMersTemp) {
//                        System.err.println(pairMerTmp.getClipLeft()+pairMerTmp.decodeCore(k-1));
//                        System.err.println(pairMerTmp.decodeCore(k-1)+pairMerTmp.getClipRight());
//                    }
                        toReport.add(terminal1.getPairMer().getPairMerString(k) + " <- terminal1");
                        toReport.add(terminal2.getPairMer().getPairMerString(k) + " <- terminal2");
                        toReport.add(pairMer.getPairMerString(k) + " <- pairMer");
//                        toReport.add(traverse(terminal1, k, true, null) + " <- traverse1");
//                        toReport.add(traverse(terminal2, k, true, null) + " <- traverse2");
                        if (previous != null) {
                            toReport.add(previous.getPairMerString(k) + " <- previous");
                        }
                        if (next != null) {
                            toReport.add(next.getPairMerString(k) + " <- next");
                        }
//                        toReport.add(pairMerNode.pairMer.getPairMerString(k) + " <-- PM String");
//                        toReport.add(pairMerNode.next.getPairMerString(k) + " <-- PM next, RC=" + pairMerNode.nextRc);
//                        toReport.add(pairMerNode.previous.getPairMerString(k) + " <-- PM previous, RC=" + pairMerNode.previousRc);
//                        toReport.add(traverse(pairMerNode, k, true, null) + " <- terminal2");
//                        toReport.add(" -=-=- ");

                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    Reporter.writeToFile(DEBUG_FILE, toReport, true);
                }
//                System.exit(1);

            }
        } else if (previous == null && next == null) {
            if (singletonNode == null) {
                singletonNode = pairMerNode;
            } else {
                Reporter.report("[BUG?]", "Another singleton node among supposedly connected PairMers?  ", getClass().getSimpleName());
            }

        }
        PairMerNode put = pairMerNodes.put(pairMer, pairMerNode);
        if (put != null) {
            Reporter.report("[BUG?]", "PairMer/Node being added agian to Map?  ", getClass().getSimpleName());
        }
    }

    /**
     * Generates the extended string assuming that connectPairMers was executed
     * earlier
     *
     * @param k
     * @return
     */
    public CharSequence toCharSeq(int k) {
        if (extendedString == null) {
            if (singletonNode != null) {
                extendedString = singletonNode.getPairMer().getPairMerString(k);
            } else if (!hasTerminalOrSingletonNode()) {
                Reporter.report("[WARNING]", "No terminal PairMerNodes ", getClass().getSimpleName() + " " + Thread.currentThread().getName());
//                for (Map.Entry<PairMer, PairMerNode> entry : pairMerNodes.entrySet()) {
//                    PairMer key = entry.getKey();
//                    PairMerNode value = entry.getValue();
//                    System.err.println(key.getPairMerString(45, "_"));
//                }               
                return null;
            } else if (terminal1 != null) {
                extendedString = traverse(terminal1, k, true, null);
            } else if (terminal2 != null) {
                extendedString = traverse(terminal2, k, true, null);
            } else {
                //not possible everything correct           
                Reporter.report("[BUG?]", "No terminal or sungleton PairMerNodes in ConnectedPairMers? ", getClass().getSimpleName());
                return null;
            }
        }
        return extendedString;
    }

    public boolean isPotentialDuplicate() {
        return potentialDuplicate;
    }

    /**
     * If a non-empty set of connected PairMers has no terminal nodes, it
     * indicates that it represents a circular molecule. In most contexts it
     * would indicate an error caused e.g. by a k-mer's reverse complement
     * accidentally matching another k-mer in the set.
     *
     * @return true if has at least one terminal node
     */
    public boolean hasTerminalOrSingletonNode() {
        return terminal1 != null || terminal2 != null || singletonNode != null;
    }

    private CharSequence traverse(PairMerNode currentNode, int k, boolean isInitial, PairMerNode alreadyVisited) {
        StringBuilder sb = new StringBuilder();
        String currentString = new String();
        if (isInitial) {
            currentString = currentNode.getPairMer().getPairMerString(k);
            isInitial = false;
        }
//        System.err.println(currentNode.getPairMer().getPairMerString(k) + " <-current");

        if (!currentNode.hasPrevious() && !currentNode.hasNext()) {
            sb.append(currentString);
        }

        if (currentNode.hasPrevious()) {
//            System.err.println(currentNode.getPrevious().getPairMerString(k) + " <-previous ");
            PairMerNode previous = pairMerNodes.get(currentNode.getPrevious());
            if (alreadyVisited == null || previous != alreadyVisited) {
                if (currentNode.isPreviousRc()) {
//                    System.err.println(SequenceOps.getReverseComplementString(previous.getPairMer().getPairMerString(k))+" <- (RC) recursion to previous");
                    sb.append(SequenceOps.getReverseComplementString(traverse(previous, k, false, currentNode)));
                    sb.append(SequenceOps.complement(currentNode.getPrevious().getClipRight()));
                } else {
//                    System.err.println(previous.getPairMer().getPairMerString(k)+" <- recursion to previous");
                    sb.append(traverse(previous, k, false, currentNode));
                    sb.append(currentNode.getPrevious().getClipLeft());
                }
                sb.append(currentString);
            } else {
            }
        }

        if (currentNode.hasNext()) {
//            System.err.println(currentNode.getNext().getPairMerString(k) + " <-next ");
            PairMerNode next = pairMerNodes.get(currentNode.getNext());
            if (alreadyVisited == null || next != alreadyVisited) {
                sb.append(currentString);
                if (currentNode.isNextRc()) {
//                    System.err.println(SequenceOps.getReverseComplementString(next.getPairMer().getPairMerString(k))+" <- (RC) recursion to next");
                    sb.append(SequenceOps.complement(currentNode.getNext().getClipLeft()));
                    sb.append(SequenceOps.getReverseComplementString(traverse(next, k, false, currentNode)));
                } else {
//                    System.err.println(next.getPairMer().getPairMerString(k)+" <- recursion to next");
                    sb.append(currentNode.getNext().getClipRight());
                    sb.append(traverse(next, k, false, currentNode));
                }
            } else {
            }
        }
//            return "not implemented yet " + currentString;

//        return sb.toString();
        return sb;
    }

    public int size() {
        return this.pairMerNodes.size();
    }

    public Set<PairMer> getKeys() {
        return pairMerNodes.keySet();
    }

//    public PairMer getTerminal1() {
//        return terminal1.getPairMer();
//    }
//
//    public PairMer getTerminal2() {
//        return terminal2.getPairMer();
//    }
    public ArrayList<PairMer> terminalMersOrSingleton() {
        ArrayList<PairMer> list = new ArrayList<>(2);
        if (terminal1 != null) {
            list.add(terminal1.getPairMer());
        }
        if (terminal2 != null) {
            list.add(terminal2.getPairMer());
        }
        if (singletonNode != null) {
            list.add(singletonNode.getPairMer());
        }
        return list;
    }

    /**
     * Class facilitates generation of a graph connecting matching PairMers
     */
    private class PairMerNode {

        PairMer pairMer;
        PairMer previous;
        boolean previousRc;
        PairMer next;
        boolean nextRc;

        public PairMerNode(PairMer pairMer, PairMer previous, boolean previousRc, PairMer next, boolean nextRc) {
            this.pairMer = pairMer;
            this.previous = previous;
            this.previousRc = previousRc;
            this.next = next;
            this.nextRc = nextRc;
        }

        public void addPrevious(PairMer previous, boolean isPreviousRC) {
            this.previous = previous;
            this.previousRc = isPreviousRC;
        }

        public void addNExt(PairMer next, boolean isNextRC) {
            this.next = next;
            this.nextRc = isNextRC;
        }

        public PairMer getPairMer() {
            return pairMer;
        }

        public PairMer getPrevious() {
            return previous;
        }

        public boolean isPreviousRc() {
            return previousRc;
        }

        public PairMer getNext() {
            return next;
        }

        public boolean isNextRc() {
            return nextRc;
        }

        public boolean hasNext() {
            return next != null;
        }

        public boolean hasPrevious() {
            return previous != null;
        }

    }

    public boolean hasAmbiguousEndTMP() {
        return (terminal1!=null && terminal1.getPairMer().isNextToAmbiguous()) || (terminal2 != null && terminal2.getPairMer().isNextToAmbiguous()) || (singletonNode !=null && singletonNode.getPairMer().isNextToAmbiguous());
    }

    public int countAmbiguousEnds() {
        int a = 0;
        if(terminal1!=null && terminal1.getPairMer().isNextToAmbiguous()) {
            a++;
        }
        if(terminal2 != null && terminal2.getPairMer().isNextToAmbiguous()) {
            a++;
        }
        if(singletonNode !=null && singletonNode.getPairMer().isNextToAmbiguous()) {
            a=9;
        }
        return a;
    }
    
}
