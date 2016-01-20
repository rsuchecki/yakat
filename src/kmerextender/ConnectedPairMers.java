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

import shared.SequenceOps;
import shared.Reporter;
import java.util.HashMap;
import java.util.Set;

/**
 * Class for storing a set of connected PairMers with the connections between them
 * 
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class ConnectedPairMers {
//    private ArrayList<PairMerNode> connectedPairMers = new ArrayList<>();

    private PairMerNode terminal1;
    private PairMerNode terminal2;
    private PairMerNode singletonNode;
    private final HashMap<PairMer, PairMerNode> pairMerNodes = new HashMap<>();

    /**
     * Given a PairMer recursively follows the implicit connections to generate a slightly more explicit 
     * graph structure, which can than be quickly traversed to generate the extended string using toString(int k)
     * @param pairMer
     * @param k
     * @param pairMersMap 
     */
    public void connectPairMers(PairMer pairMer, int k, PairMersMap pairMersMap) {
        if (!pairMer.isVisited()) {
            pairMer.setVisited();
            PairMerStrings wrapper = new PairMerStrings(pairMer, k);
            //
            String otherCoreOfKmer1 = wrapper.getOtherCoreOfKmer1();
            PairMer otherPairMer1 =  pairMersMap.get(otherCoreOfKmer1, k);
            String otherCoreOfKmer2 = wrapper.getOtherCoreOfKmer2();
            PairMer otherPairMer2 =  pairMersMap.get(otherCoreOfKmer2, k);
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
                if (!otherCoreOfKmer1.equals(otherPairMer1.decodeCore(k - 1))) {
                    otherPairMer1isRC = true;
                }
                if (!otherPairMer1.isVisited()) {
                    connectPairMers(otherPairMer1, k, pairMersMap);
//                    terminalPairMer1 = addPairMerToListOfConnected(otherPairMer1, k, pairMersMap, connectedPairMers);
                }
            }
            boolean otherPairMer2isRC = false;
//            System.err.println(wrapper.getKmer2String() + "\tk2");
//            System.err.println(SequenceOps.getReverseComplementString(wrapper.getKmer2String()) + "\tk2 RC");
            if (otherPairMer2 != null) {
//                System.err.println(wrapper.getOtherCoreOfKmer2() + "\tk2otherCore");
//                System.err.println(SequenceOps.getReverseComplementString(wrapper.getOtherCoreOfKmer2()) + "\tk2otherCoreRC");
                if (!otherCoreOfKmer2.equals(otherPairMer2.decodeCore(k - 1))) {
                    otherPairMer2isRC = true;
                }
                if (!otherPairMer2.isVisited()) {
                    connectPairMers(otherPairMer2, k, pairMersMap);
                }
            }
            add(pairMer, otherPairMer1, otherPairMer1isRC, otherPairMer2, otherPairMer2isRC);

        }
    }
    
    /**
     * 
     * @param pairMer
     * @param previous
     * @param previousRc
     * @param next
     * @param nextRc 
     */
    private void add(PairMer pairMer, PairMer previous, boolean previousRc, PairMer next, boolean nextRc) {
        PairMerNode pairMerNode = new PairMerNode(pairMer, previous, previousRc, next, nextRc);
//        connectedPairMers.add(pairMerNode);
        if (previous == null ^ next == null) {
            if (terminal1 == null) {
                terminal1 = pairMerNode;
            } else if (terminal2 == null) {
                terminal2 = pairMerNode;
            } else {                
                Reporter.report("[BUG?]", "Third terminal PairMerNode identified  ", getClass().getSimpleName() );
//                System.err.println(terminal1.getPairMer().getPairMerString(30)+" <- terminal1");
//                System.err.println(terminal2.getPairMer().getPairMerString(30)+" <- terminal2");
//                System.err.println(pairMer.getPairMerString(30)+" <- pairMer");
////                System.err.println(previous.getPairMerString(30)+" <- previous");
////                System.err.println(next.getPairMerString(30)+" <- next");
//                System.err.println(" -=-=- ");
                
            }
        } else if (previous == null && next == null) {
            if (singletonNode == null) {
                singletonNode = pairMerNode;
            } else {
                Reporter.report("[BUG?]", "Another singleton node among supposedly connected PairMers?  ", getClass().getSimpleName() );
            }

        }
        PairMerNode put = pairMerNodes.put(pairMer, pairMerNode);
        if (put != null) {
            Reporter.report("[BUG?]", "PairMer/Node being added agian to Map?  ", getClass().getSimpleName() );
        }
    }

    /**
     * Generates the extended string assuming that connectPairMers was executed earlier
     * @param k
     * @return 
     */
    public String toString(int k) {
        if (singletonNode != null) {
            return singletonNode.getPairMer().getPairMerString(k);
        } else if (!hasTerminalOrSingletonNode()) { 
            Reporter.report("[WARNING]", "No terminal PairMerNodes ", getClass().getSimpleName() );
            return null;
        } else if (terminal1 != null) {
            return traverse(terminal1, k, true, null);
        } else if (terminal2 != null) {
            return traverse(terminal2, k, true, null);
        } else {
            //not possible            
            Reporter.report("[BUG?]", "No terminal or sungleton PairMerNodes in ConnectedPairMers? ", getClass().getSimpleName() );
            return null;
        }

    }
    
    /**
     * If a non-empty set of connected PairMers has no terminal nodes, it 
     * indicates that it represents a circular molecule. In most contexts 
     * it would indicate an error caused e.g. by a k-mer's reverse complement 
     * accidentally matching another k-mer in the set. 
     * @return true if has at least one terminal node 
     */
    public boolean hasTerminalOrSingletonNode() {
        return terminal1 != null || terminal2 != null || singletonNode != null;
    }
 
    private String traverse(PairMerNode currentNode, int k, boolean isInitial, PairMerNode alreadyVisited) {
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

        return sb.toString();
    }

    public int size() {
        return this.pairMerNodes.size();
    }

    public Set<PairMer> getKeys() {
        return pairMerNodes.keySet();
    }

    public PairMer getTerminal1() {
        return terminal1.getPairMer();
    }

    public PairMer getTerminal2() {
        return terminal2.getPairMer();
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
}
