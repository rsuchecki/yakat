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

//import gnu.trove..set.hash.THashSet;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ConcurrentNavigableMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.atomic.AtomicLong;
import shared.Reporter;

/**
 * Wrapper around a concurrent collection (ConcurrentSkipListMap, but other
 * options possible) of PairMer objects. After the Map is populated, use
 * .purge()
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class PairMersMap extends shared.MerMap {

//    private final Object LOCK;
    private ConcurrentSkipListMap<PairMer, PairMer> pairMersSkipListMap;
    private ConcurrentSkipListMap<PairMer, PairMer> terminalPairMers;
//    private ConcurrentSkipListMap<PairMer, PairMer> pairMersSkipListMap;
    private Integer k;
    private PairMersMap parentMap;
    private AtomicLong ambiguous;
    private AtomicLong size;
    private AtomicLong sizeTerminal = new AtomicLong();

////    private PairMer[] prefixToPairMers;
//    private ChronicleMap<long[], Long> chronicleMap;
//    BTreeMap<long[], Long> mapDB;

    /**
     * Instantiate the Map
     *
     * @param k
     */
    public PairMersMap(Integer k) {
//        prefixToPairMers = new PairMer[16777216];
        this.size = new AtomicLong();
        this.ambiguous = new AtomicLong();
        this.sizeTerminal = new AtomicLong();
        if (k != null) {
            pairMersSkipListMap = new ConcurrentSkipListMap<>();
            terminalPairMers = new ConcurrentSkipListMap<>();
//           pairMersHashMap  = new ConcurrentHashMap<>();
            this.k = k;
        }
    }

//    public void initChronicleMap() {
//        System.err.println("Hardcoded tests of chronicle map @ k=71 ...");
//        long[] encodedCoreLongArray = CoreCoder.encodeCoreLongArray("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
//        ChronicleMapBuilder<long[], Long> chronicleMapBuilder
//                = ChronicleMapBuilder.of(long[].class, Long.class)
//                        .name("chronicleMap")
//                        .constantKeySizeBySample(encodedCoreLongArray)
//                        .entries(5_000_000_000L);
////                        .entries(2_000_000L);
//        chronicleMap = chronicleMapBuilder.create();
//
//        System.err.println("Hardcoded tests of chronicle map @ k=71 - map created");
////        chronicleMap.put(encodedCoreLongArray, Long.MIN_VALUE);
////        System.err.println("Hardcoded tests of chronicle map @ k=71 - test element added");
//    }
//
//    public void initMapDB() {
//        System.err.println("Hardcoded tests of mapDB  ...");
//        DB db = DBMaker
//            .memoryDirectDB()
////                .fileDB("file.db")
////                .fileMmapEnable()
////                .allocateStartSize(10_000_000)
//                .transactionEnable()
//                .make();
//        mapDB = db.treeMap("map",Serializer.LONG_ARRAY, Serializer.LONG)
//                .createOrOpen();
//        System.err.println("Hardcoded tests of mapDB  ... created");
//        
//        
//        // using native order speeds access for values longer than a byte.
////ByteBuffer bb = ByteBuffer.allocateDirect(1024*1024*1024).order(ByteOrder.nativeOrder());
////// start at some location.
////bb.position(0);
////bb.put((byte) 1);
////bb.putInt(6456456);
////bb.putDouble(1.4564687);
////
////// to read back.
////bb.position(0);
////byte b = bb.get();
////int i = bb.getInt();
////double d = bb.getDouble();
////        
////int x =0;
////        BTreeMap<byte[], Long> map = db.treeMap("map")
////                .keySerializer(Serializer.BYTE_ARRAY)
////                .valueSerializer(Serializer.LONG)
////                .createOrOpen();
//    }
//    
//    public void initOHC() {
//        OHCache ohCache = OHCacheBuilder.newBuilder()
//                                .keySerializer(new CacheSerializer<Object>() {
//            @Override
//            public void serialize(Object t, ByteBuffer bb) {
//                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//            }
//
//            @Override
//            public Object deserialize(ByteBuffer bb) {
//                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//            }
//
//            @Override
//            public int serializedSize(Object t) {
//                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//            }
//        })
//                                .valueSerializer(Serializer.LONG)
//                                .build();
//    }

    


    /**
     * Constructor for exposing a subMap for multithreaded purging
     *
     * @param k
     * @param subMap
     * @param parentMap
     */
    public PairMersMap(Integer k, ConcurrentNavigableMap<PairMer, PairMer> subMap, PairMersMap parentMap) {
        this.k = k;
        this.pairMersSkipListMap = new ConcurrentSkipListMap<>(subMap);
        this.terminalPairMers = parentMap.getTerminalPairMers();
        this.sizeTerminal = parentMap.getSizeTerminal();
        this.ambiguous = parentMap.getAmbiguous();
        this.size = parentMap.getSize();
        this.parentMap = parentMap;
    }

    public AtomicLong getAmbiguous() {
        return ambiguous;
    }

    public AtomicLong getSize() {
        return size;
    }

    public AtomicLong getSizeTerminal() {
        return sizeTerminal;
    }

    public boolean isEmpty() {
        return pairMersSkipListMap == null || pairMersSkipListMap.isEmpty();
    }

    public boolean isNull() {
        return pairMersSkipListMap == null;
    }

//    /**
//     * First tries to atomically add a k-mer to the Map, if this fails,
//     * synchronized method is used to update the previously stored PairMer
//     *
//     * @param kmerString
//     * @param frontClip
//     * @param overlapLength : normally k-1
//     * @param inputKmersUnique : true when a list of unique k-mers given as
//     * input, false if k-mers extracted directly from FASTA/FASTQ
//     */
//    public void addToPairMersMap(String kmerString, boolean frontClip, int overlapLength, boolean inputKmersUnique) {
//        PairMer pairMer = PairMerGenerator.generatePairMer(kmerString, frontClip, overlapLength);
//
//        //Atomic operation START
//        PairMer previousStoredPairMer = pairMersSkipListMap.putIfAbsent(pairMer, pairMer);
//        //Atomic operation END
//        if (previousStoredPairMer != null) {
//            //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
//            previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique);
//        }
//    }
//    Random r = new Random();

    /**
     * First tries to atomically add a k-mer to the Map, if this fails,
     * synchronized method is used to update the previously stored PairMer
     *
     * @param sequence
     * @param from
     * @param to
     * @param frontClip
     * @param inputKmersUnique
     * @param freq
     */
    public void addToPairMersMap(CharSequence sequence, int from, int to, boolean frontClip, boolean inputKmersUnique, int freq) {
//        synchronized (KmerExtender.class) {
        boolean failed = false;
        PairMer pairMer = null;
        try {
            pairMer = PairMerTypeSelector.generatePairMer(sequence, from, to, frontClip, freq);

//            String canonicalKmer;
//            if(pairMer.hasLeftClip()) {
//                canonicalKmer = SequenceOps.getCanonical(pairMer.getClipLeft()+pairMer.decodeCore(44));                
//            } else {
//                canonicalKmer = SequenceOps.getCanonical(pairMer.decodeCore(44)+pairMer.getClipRight());                                
//            }
//            String canonicalDecoded = SequenceOps.getCanonical(sequence.subSequence(from, to+1).toString());
//            if(!canonicalDecoded.equals(canonicalKmer)) {
//                System.err.println("Kmer - decode mismtach "+canonicalKmer+" "+canonicalDecoded);
//            }
        } catch (NonACGTException ex) {
//            Reporter.report("[WARNING]", ex.getMessage(), getClass().getCanonicalName());
            failed = true;
        }
        if (!failed) {
//            ConcurrentSkipListSet<PairMer> pairMersSkipListSet = new ConcurrentSkipListSet<>();
//            pairMersSkipListSet.

            //Atomic operation START                    
            PairMer previousStoredPairMer = pairMersSkipListMap.putIfAbsent(pairMer, pairMer);
            //Atomic operation END
//            System.err.println("Putting elem...");
//            chronicleMap.putIfAbsent(((PairMer3LongEncoded) pairMer).getBitFields(), r.nextLong());
//            mapDB.putIfAbsent(((PairMer2LongEncoded) pairMer).getBitFields(), r.nextLong());
//            System.err.println("Done");
//            chronicleMap.

            if (previousStoredPairMer == null) {
                size.incrementAndGet();
            } else {
                //TODO!!!! If rc of a seq == seq we don't wan't duplciates i.e. 2 clipmers derived from a single kmer[checking the underlying kmer not just the clipped part]
                previousStoredPairMer.addKmerSynchronized(pairMer, inputKmersUnique, freq);
            }
        }
//        }

    }

    /**
     * Given a core string, retrieves the matching PairMer
     *
     * @param core, which will be converted to its canonical form
     * @param k
     * @return PairMer if present in Map, null otherwise
     * @throws kmerextender.NonACGTException
     */
    public PairMer get(CharSequence core, int k) throws NonACGTException {
        return pairMersSkipListMap.get(PairMerTypeSelector.getPairMer(core, k));
    }

//    /**
//     * Given a core string, retrieves the matching PairMer from the list of
//     * identified terminal PMs
//     *
//     * @param core, which will be converted to its canonical form
//     * @param k
//     * @return PairMer if present in Map, null otherwise
//     * @throws kmerextender.NonACGTException
//     */
//    public PairMer getTerminal(CharSequence core, int k) throws NonACGTException {
//        return terminalPairMers.get(PairMerGenerator.getPairMer(core, k));
//    }
    /**
     * Retrieve a PairMer with a matching core
     *
     * @param anotherPairMer
     * @return
     */
    public PairMer get(PairMer anotherPairMer) {
        return pairMersSkipListMap.get(anotherPairMer);
    }

    public long size() {
        return size.longValue();
    }

    public long sizeTerminals() {
        return sizeTerminal.longValue();
    }

    public boolean remove(PairMer elem) {
        if (pairMersSkipListMap.remove(elem) != null) {
            size.decrementAndGet();
            return true;
        }
        return false;
    }

    /**
     *
     * @return the underlying data structure
     */
    public ConcurrentSkipListMap getPairMersMap() {
        return pairMersSkipListMap;
    }

    public boolean contains(PairMer key) {
        return pairMersSkipListMap.containsKey(key);
    }

    public Iterator<PairMer> iterator() {
        return pairMersSkipListMap.keySet().iterator();
    }

    public ConcurrentSkipListMap<PairMer, PairMer> getTerminalPairMers() {
        return terminalPairMers;
    }

//    private boolean consistentTraversal(PairMer current) {
//        PairMer left = pairMersSkipListMap.get(current.getOtherPairmerCoreLeft(k));
//        PairMer right = pairMersSkipListMap.get(current.getOtherPairmerCoreRight(k));
//        PairMer otherPairmerCoreLeftLeft = left.getOtherPairmerCoreLeft(k);
//        PairMer otherPairmerCoreLeftRight = left.getOtherPairmerCoreRight(k);
//        PairMer otherPairmerCoreRightLeft = right.getOtherPairmerCoreLeft(k);
//        PairMer otherPairmerCoreRightRight = right.getOtherPairmerCoreRight(k);
//        if (!current.equals(pairMersSkipListMap.get(otherPairmerCoreLeftLeft)) && !current.equals(pairMersSkipListMap.get(otherPairmerCoreLeftRight))
//                && !current.equals(pairMersSkipListMap.get(otherPairmerCoreRightLeft)) && !current.equals(pairMersSkipListMap.get(otherPairmerCoreRightRight))) {
//                        return false;
//        }
//        return true;
//    }
//
//    private void printPairMerDetails(String prefix, PairMer pm) {
//        System.err.println(prefix+"\t"+ pm.getPairMerString(k) + "\t" + pm.getStoredCountLeft() + "\t" + pm.getStoredCountRigth() + "\t" + pm.isInvalid() + "\tPURGED");
//
//    }
    /**
     * Removes from the Map each PairMer which (i) represents ambiguous
     * extension (>2 k-mers matching the core) or (ii) has no extensions (only
     * one k-mer matching the core) or (iii) represents 2 or more conflicting
     * k-mers (matching core, but both extending in the same direction) Run only
     * after the set/map is fully populated.
     *
     * Could be parallelized by applying the purge method to subsets see
     * .subSet()
     *
     * @param minKmerFrequency
     * @param threads
     * @return the number of elements removed
     */
    public long purge(int minKmerFrequency) {

        long count = 0L;

        Iterator<PairMer> it = pairMersSkipListMap.keySet().iterator();

//        Reporter.report("[WARNING]", "Pointless map traversal started", getClass().getCanonicalName());
//        long c = 0;
//        while (it.hasNext()) {
//            PairMer next = it.next();
//            c++;
//        }
//        Reporter.report("[WARNING]", "Pointless map traversal finished, n="+NumberFormat.getInstance().format(c), getClass().getCanonicalName());
//        it = pairMersSkipListMap.keySet().iterator();
//        PairMer firstKey = pairMersSkipListMap.firstKey();
//        ArrayList<ConcurrentNavigableMap<PairMer, PairMer>> mapChunks = new ArrayList<>();
//        ArrayList<PairMer> spacedPairMers = firstKey.generateSpacedPairMersForMapSplitting(k, 10);
//        Random r = new Random();
//        System.err.println("Running recursive map splitting");
//        if (threads > 1) {
//            long attempts = recursiveSplitMap(pairMersSkipListMap, mapChunks, size() / threads / 100, size() / threads / 10);
//            System.err.println("Finished running recursive map splitting, total attempts = " + attempts);
//        }
//        int[] sizes = new int[mapChunks.size()];
//        long total = 0;
//        long min = Long.MAX_VALUE;
//        long max = Long.MIN_VALUE;
//        for (int j = 0; j < mapChunks.size(); j++) {
////            System.err.println("chunk " + j + " size = " + mapChunks.get(j).size());
//            sizes[j] = mapChunks.get(j).size();
//            total += mapChunks.get(j).size();
//            min = sizes[j] < min ? sizes[j] : min;
//            max = sizes[j] > max ? sizes[j] : max;
//        }
//        System.err.println(sizes.length + " chunks, tot=" + NumberFormat.getNumberInstance().format(total) + ", mean=" + total / sizes.length + ", median=" + sizes[sizes.length / 2] + ", min=" + min + ", max=" + max);
        while (it.hasNext()) {
            PairMer next = it.next();
//            System.err.println(next.toString());

//            PairMerStrings wrapper = new PairMerStrings(next, k);
//            //
//            String otherCoreOfKmer1 = wrapper.getOtherCoreOfKmer1();
//            PairMer otherPairMer1 = null;
//            String otherCoreOfKmer2 = wrapper.getOtherCoreOfKmer2();
//            PairMer otherPairMer2 = null;
//            try {
//                otherPairMer1 = get(otherCoreOfKmer1, k);
//                PairMerStrings wrapper1 = new PairMerStrings(otherPairMer1, k);
//                String otherCoreOfKmer11 = wrapper1.getOtherCoreOfKmer1();
//                PairMer otherPairMer11 = get(otherCoreOfKmer11, k);
//                String otherCoreOfKmer12 = wrapper1.getOtherCoreOfKmer2();
//                PairMer otherPairMer12 = get(otherCoreOfKmer12, k);
//                
//                
//                otherPairMer2 = get(otherCoreOfKmer2, k);
//                PairMerStrings wrapper2 = new PairMerStrings(otherPairMer2, k);
//                String otherCoreOfKmer21 = wrapper2.getOtherCoreOfKmer1();
//                PairMer otherPairMer21 = get(otherCoreOfKmer21, k);
//                String otherCoreOfKmer22 = wrapper2.getOtherCoreOfKmer2();
//                PairMer otherPairMer22 = get(otherCoreOfKmer22, k);
//                
//                if(!next.equals(otherPairMer11) && !next.equals(otherPairMer12)
//                    && !next.equals(otherPairMer21) && !next.equals(otherPairMer22)) {
//                    Reporter.report("[BUG]", "You'd expect that one of the pairmers linked to a pairmer would be that pairmer", getClass().getCanonicalName());
//                    
//                }
//                
//            } catch (NonACGTException ex) {
//                Reporter.report("[WARNING]", "Unexpected NonACGTException caught", getClass().getCanonicalName());
//            }
//            String decodeCore = next.decodeCore(44);
//            System.err.println(next.getPairMerString(k)+"\t"+next.getStoredCountLeft()+"\t"+next.getStoredCountRigth()+"\t"+next.isInvalid());
//            if (next.isInvalid()) {
//                it.remove();
//                count++;
            if (next.isInvalid() || next.getStoredCountLeft() == 0 || next.getStoredCountRigth() == 0) {
//                //PM removed either due to ambig or holding just one k-mer, any valid adjacent k-mer present in map must be a terminal one
                if (next.isInvalid()) {
                    ambiguous.incrementAndGet();
                }
                try {
                    String decodedCore = next.decodeCore(k - 1);
                    if (next.hasLeftClip()) {
                        StringBuilder otherCoreOfKmer1 = new StringBuilder();
                        otherCoreOfKmer1.append(next.getClipLeft());
                        otherCoreOfKmer1.append(decodedCore.subSequence(0, decodedCore.length() - 1));
                        if (addTerminal(otherCoreOfKmer1, minKmerFrequency)) {
                            sizeTerminal.incrementAndGet();
                        }

                    }
                    if (next.hasRightClip()) {
                        StringBuilder otherCoreOfKmer2 = new StringBuilder();
                        otherCoreOfKmer2.append(decodedCore.subSequence(1, decodedCore.length()));
                        otherCoreOfKmer2.append(next.getClipRight());
                        if (addTerminal(otherCoreOfKmer2, minKmerFrequency)) {
                            sizeTerminal.incrementAndGet();
                        }
                    }
                } catch (NonACGTException ex) {
                    Reporter.report("[WARNING]", "Unexpected NonACGTException caught", getClass().getCanonicalName());
                }
                getParentMap().remove(next);
                size.decrementAndGet();
                count++;
//                }
//            } else if (!terminalPairMers.containsKey(next)) { //if not already marked as terminal
//                try {
//                    String decodedCore = next.decodeCore(k - 1);
//                    if (next.hasLeftClip()) {
//                        StringBuilder otherCoreOfKmer1 = new StringBuilder();
//                        otherCoreOfKmer1.append(next.getClipLeft());
//                        otherCoreOfKmer1.append(decodedCore.subSequence(0, decodedCore.length() - 1));
//                        PairMer merInMap = get(PairMerGenerator.getPairMer(otherCoreOfKmer1, k));
//                        if(merInMap == null) {
//                            terminalPairMers.putIfAbsent(next, next);
//                            continue; //so that we dont compute the other core just to place it among terminal pair mers again if is singleton
//                        }
//
//                    }
//                    if (next.hasRightClip()) {
//                        StringBuilder otherCoreOfKmer2 = new StringBuilder();
//                        otherCoreOfKmer2.append(decodedCore.subSequence(1, decodedCore.length()));
//                        otherCoreOfKmer2.append(next.getClipRight());
//                        PairMer merInMap = get(PairMerGenerator.getPairMer(otherCoreOfKmer2, k));
//                        if(merInMap == null) {
//                            terminalPairMers.putIfAbsent(next, next);
//                        }                        
//                    }
//                } catch (NonACGTException ex) {
//                    Reporter.report("[WARNING]", "Unexpected NonACGTException caught", getClass().getCanonicalName());
//                }
            }
//else if (next.getStoredCountLeft() < minKmerFrequency || next.getStoredCountRigth() < minKmerFrequency) {
////                System.err.println("\tpurge\tisInvalid="+next.isInvalid());
////            if (next.isInvalid() || next.getStoredCount() != 2) {
////            if (next.isInvalid() || next.getStoredCount() < 2) {
////                System.err.println("Removing:");//\n\t"+next.getClipLeft()+"-"+next.getTmpCore()+"-"+next.getClipRight());
////                for (String s : next.getHistory()) {
////                    System.err.println("\t" + s);
////                }
//
//                it.remove();
////                pairMersSkipListMap.remove(next);
//                count++;
////                if (next.hasBothClips()) {
////                    System.err.println(next.getPairMerString(k) + "\t" + next.getStoredCountLeft() + "\t" + next.getStoredCountRigth() + "\t" + next.isInvalid() + "\tPURGED");
//////                    System.err.println(next.getPairMerString(k) + " <- PURGED");
////                }
////            } else {
////                System.err.println("\tkeep");
////                System.err.println(next.getPairMerString(45)+"\tREMOVED");
//
////            } else {
////                System.err.println(next.getPairMerString(45)+"\tKEPT");
//            }
        }

//        //EXPERIMENTS only
//        it = pairMersSkipListMap.keySet().iterator();
//        int buckets = pairMersSkipListMap.size();
//        byte[] hashKeys = new byte[buckets];
//        while (it.hasNext()) {
//            PairMer next = it.next();
//            int hashcode = next.hashCode() % buckets;
//            if (hashcode < 0) {
//                hashKeys[-hashcode]++;
//            } else {
//                hashKeys[hashcode]++;
//            }
//            total++;
//        }
//        double mean = (double) total / buckets;
//        double deviation = 0;
//        int single = 0;
//        int multi = 0;
//        int empty = 0;
//        for (int i = 0; i < hashKeys.length; i++) {
////            System.err.printf("%8d %8d\n",i,hashKeys[i]);
//            deviation += Math.pow((hashKeys[i] - mean), 2);
//            if (hashKeys[i] < 1) {
//                empty++;
//            } else if (hashKeys[i] == 1) {
//                single++;
//            } else {
//                multi++;
//            }
//        }
//        System.err.println("Mean  = " + mean);
//        System.err.println("StDev = " + Math.sqrt(deviation / buckets));
//        System.err.println("Empty = " + empty + " " + ((double) empty / buckets) + "%");
//        System.err.println("Single= " + single + " " + ((double) single / buckets) + "%");
//        System.err.println("Multi = " + multi + " " + ((double) multi / buckets) + "%");
////        for (int i = hashKeys.length-10000; i < hashKeys.length; i++) {
////            System.err.printf("%8d %8d\n",i,hashKeys[i]);
////            
////        }
        return count;
    }

    private boolean addTerminal(CharSequence pairMerCore, int minKmerFrequency) throws NonACGTException {
        //                        PairMer encodedCoreOfKmer2 = PairMerGenerator.getPairMer(otherCoreOfKmer2, k);
////                        System.err.println(otherCoreOfKmer2.toString() + "\t CREATED FROM DECODED CORE (2)");
////                synchronized (KmerExtender.class) {
////                        System.err.println(next.getPairMerString(k) + "\tREMOVING");
//
//                        PairMer otherPairMer2 = getParentMap().get(encodedCoreOfKmer2);
//                        if (otherPairMer2 != null && !otherPairMer2.isInvalid() && otherPairMer2.getStoredCountLeft() >= minKmerFrequency && otherPairMer2.getStoredCountRigth() >= minKmerFrequency) {
////                        pairMersSkipListMap.remove(otherPairMer2);
////                            System.err.println(otherPairMer2.getPairMerString(k) + "\tADDING TO TERMINAL");
//                            PairMer put = terminalPairMers.putIfAbsent(otherPairMer2, otherPairMer2);
//                            if (put != null) {
////                                System.err.println(otherPairMer2.getPairMerString(k) + "\t - already in");
//                                int x = 0;
//                            }
//                        }
        PairMer encodedCoreOfKmer1 = PairMerTypeSelector.getPairMer(pairMerCore, k);
        PairMer otherPairMer1 = getParentMap().get(encodedCoreOfKmer1);
        if (otherPairMer1 != null && !otherPairMer1.isInvalid() && otherPairMer1.getStoredCountLeft() >= minKmerFrequency && otherPairMer1.getStoredCountRigth() >= minKmerFrequency) {
            //                        pairMersSkipListMap.remove(otherPairMer1);
//                            System.err.println(otherPairMer1.getPairMerString(k) + "\tADDING TO TERMINAL");
            return terminalPairMers.putIfAbsent(otherPairMer1, otherPairMer1) == null;
        }
        return false;
    }

    //    public boolean isOutOfMemory() {
    //        return OutOfMemory;
    //    }
    //
    //    public synchronized void setOutOfMemory() {
    //        this.OutOfMemory = true;
    //    }
    //    public int purgeSeedMers(PairMersMap seedMersMap) {
    //        Iterator<PairMer> it = seedMersMap.getPairMersSkipListMap().keySet().iterator();
    //        
    //    }
    public Integer getK() {
        return k;
    }

    public PairMersMap getParentMap() {
        return parentMap == null ? this : parentMap;
    }

//    public long recursiveSplitMap(ArrayList<ConcurrentNavigableMap<PairMer, PairMer>> submaps, int minChunk, int maxChunk, ArrayBlockingQueue<PairMersMap> mapsQueue) {
//        return recursiveSplitMap(pairMersSkipListMap, submaps, minChunk, maxChunk, k, mapsQueue, this);
//    }
//
//    private long recursiveSplitMap(ConcurrentNavigableMap<PairMer, PairMer> map, ArrayList<ConcurrentNavigableMap<PairMer, PairMer>> submaps, int minChunk, int maxChunk, int k,
//            ArrayBlockingQueue<PairMersMap> mapsQueue, PairMersMap masterMap) {
//        PairMer higherKey = null;
//        long count = 0;
//        ConcurrentNavigableMap<PairMer, PairMer> headMap;
//        ConcurrentNavigableMap<PairMer, PairMer> tailMap;
//        do {
////            do {
//                try {
//                    PairMer luckyDipPairMer = map.firstKey().luckyDipPairMer(k, map.lastKey());
//                    count++;
////                System.err.println(luckyDipPairMer.toString()+ "<- lucky dip");                
//                    higherKey = map.higherKey(luckyDipPairMer);
//                    if (higherKey == null) {
//                        int x = 0;
////                    System.err.println(higherKey.toString()+ " <- higher");
//                    }
//                } catch (NullPointerException e) {
//
//                }
//            } while (higherKey == null);
//            System.err.println("Tried " + NumberFormat.getNumberInstance().format(count));
//
//            headMap = map.headMap(higherKey, false);
//            tailMap = map.tailMap(higherKey, true);
////        } while (headMap.size() < minChunk || tailMap.size() < minChunk);
////        System.err.println("Tried " + count + " keys, |HeadMap|=" + headMap.size() + " |TailMap|=" + tailMap.size() + " |submaps|=" + submaps.size());
//
//        if (headMap.size() <= maxChunk) {
//            try {
//                mapsQueue.put(new PairMersMap(k, headMap, masterMap));
//            } catch (InterruptedException ex) {
//                ex.printStackTrace();
//            }
//            submaps.add(headMap);
//        } else {
//            count += recursiveSplitMap(headMap, submaps, minChunk, maxChunk, k, mapsQueue, masterMap);
//        }
//
//        if (tailMap.size() <= maxChunk) {
//            try {
//                mapsQueue.put(new PairMersMap(k, tailMap, masterMap));
//            } catch (InterruptedException ex) {
//                ex.printStackTrace();
//            }
//            submaps.add(tailMap);
//        } else {
//            count += recursiveSplitMap(tailMap, submaps, minChunk, maxChunk, k, mapsQueue, masterMap);
//        }
//        return count;
//    }
}
