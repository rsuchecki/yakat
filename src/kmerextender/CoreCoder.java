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

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import argparser.PositionalOpt;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.logging.Level;
import java.util.logging.Logger;
import shared.InputReaderProducer;
import shared.Reporter;
import shared.SequenceOps;

/**
 * CAREFUL! - using all bits in signed fields, may cause errors if used for
 * comparisons without decoding (for example lex order of core vs its RC)
 *
 *
 *
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class CoreCoder {

    public static int[] encodeCoreIntArray(CharSequence kmerSequence) {
        int INT_LENGTH = 32; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = kmerSequence.length();
        int intsNeeded = (int) Math.ceil((double) stringLength * 2 / INT_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = intsNeeded * INT_LENGTH - stringLength * 2;
        int kmerCoreBitsArray[] = new int[intsNeeded];
//        char[] kmerCharArray = kmerSequence.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < INT_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
//                switch (kmerCharArray[position]) {
                switch (kmerSequence.charAt(position)) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    default:
                        System.err.println("Failed ecoding kmerstring to int array....");
                        System.err.println("Offending char: " + kmerSequence.charAt(position));
                        System.err.println("in " + kmerSequence);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }
        return kmerCoreBitsArray;
    }

    public static String decodeCore(int encodedSequenceLength, int kmerCoreBitsArray[]) {
        int INT_LENGTH = 32; //32 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % INT_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = INT_LENGTH - 1;
            if (i == 0 && lastChunkLength != 0) {
                startPrintingBitsFrom = lastChunkLength - 1;
            }
            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
                int b1 = kmerCoreBitsArray[i] >> j & 1;
                int b2 = kmerCoreBitsArray[i] >> j - 1 & 1;
//                int b = (kmerCoreBitsArray[i] & (1 << j)) >> j;
//                int m = (kmerCoreBitsArray[i] & (1 << j - 1)) >> j - 1;
                if (b1 == 0 && b2 == 0) {
                    sb.append("A");
                } else if (b1 == 0 && b2 == 1) {
                    sb.append("C");
                } else if (b1 == 1 && b2 == 0) {
                    sb.append("G");
                } else if (b1 == 1 && b2 == 1) {
                    sb.append("T");
                }
            }
        }
        return sb.toString();
    }

    public static long[] encodeCoreLongArray(CharSequence coreSequence) {
        int LONG_LENGTH = 64; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        int stringLength = coreSequence.length();
        int longsNeeded = (int) Math.ceil((double) stringLength * 2 / LONG_LENGTH); //number of ints needed to store this String 
//        System.out.println("Need "+intsNeeded+" "+ INT_LENGTH +" bit word(s) to encode "+kmerString.length()+ " nucl sequence");
        int currentInt = 0;
        int position = 0;
        int positionInInt = longsNeeded * LONG_LENGTH - stringLength * 2;
        long kmerCoreBitsArray[] = new long[longsNeeded];
//        char[] kmerCharArray = coreSeq.toCharArray();
        while (position < stringLength) {
            while ((positionInInt < LONG_LENGTH) && (position < stringLength)) {
                kmerCoreBitsArray[currentInt] <<= 1;
                switch (coreSequence.charAt(position)) {
                    case 'A':
                    case 'a':
                        //if A : 00
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'C':
                    case 'c':
                        //if C : 01
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    case 'G':
                    case 'g':
                        //if G : 10
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        break;
                    case 'T':
                    case 't':
                        //if T : 11
                        kmerCoreBitsArray[currentInt]++;
                        kmerCoreBitsArray[currentInt] <<= 1;
                        kmerCoreBitsArray[currentInt]++;
                        break;
                    default:
                        System.err.println("Failed ecoding kmerstring to int array....");
                        System.err.println("Offending char: " + coreSequence.charAt(position));
                        System.err.println("in " + coreSequence);
                        System.err.println("....exiting");
                        System.exit(1);
                }
                positionInInt += 2;
                position++;
            }
            currentInt++;
            positionInInt = 0;
        }
//        String decodeCore = decodeCore(stringLength, kmerCoreBitsArray);
//        if (!decodeCore.equals(coreString)) {
////            System.err.println("Error encoding/decoding " + kmerCoreBits1 + " " + kmerCoreBits2);
//            System.err.println(coreString + " <-core");
//            System.err.println(decodeCore + " <-decoded");
//            decodeCore = decodeCore(stringLength, kmerCoreBitsArray);
//        }
        return kmerCoreBitsArray;
    }

    public static String decodeCore(int encodedSequenceLength, long kmerCoreBitsArray[]) {
        int LONG_LENGTH = 64; //64 - sign bit -1 to make even as 2 bits stored per nucleotide
        StringBuilder sb = new StringBuilder();
        int lastChunkLength = encodedSequenceLength * 2 % LONG_LENGTH;
        for (int i = 0; i < kmerCoreBitsArray.length; i++) {
            int startPrintingBitsFrom = LONG_LENGTH - 1;
            if (i == 0 && lastChunkLength != 0) {
                startPrintingBitsFrom = lastChunkLength - 1;
            }
            for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
                long b1 = kmerCoreBitsArray[i] >> j & 1;
                long b2 = kmerCoreBitsArray[i] >> j - 1 & 1;
                if (b1 == 0 && b2 == 0) {
                    sb.append("A");
                } else if (b1 == 0 && b2 == 1) {
                    sb.append("C");
                } else if (b1 == 1 && b2 == 0) {
                    sb.append("G");
                } else if (b1 == 1 && b2 == 1) {
                    sb.append("T");
                }
            }
        }
        return sb.toString();
    }

//    public static int compareCores(int[] core, int[] anotherCore) {
//        for (int i = 0; i < core.length; i++) {
//            try {
//                if (core[i] < anotherCore[i]) {
//                    return -1;
//                } else if (core[i] > anotherCore[i]) {
//                    return 1;
//                }
//            } catch (IndexOutOfBoundsException e) {
//                System.err.println(e.getMessage());
//            }
//        }
//        return 0;
//    }

    public static int compareCores(long[] core, long[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
//                if (core[i] < anotherCore[i]) {
//                (longA < longB) ^ (longA < 0) ^ (longB< 0) ? 1 : -1;
                //tryuing unisgned comparison to use all bits
                if (core[i] == anotherCore[i]) {
                    //keep going
                } else if ((core[i] < anotherCore[i]) ^ (core[i] < 0) ^ (anotherCore[i] < 0)) {
                    return 1;
                } else {
                    return -1;
                }
            } catch (IndexOutOfBoundsException e) {
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }
    
    public static int compareCores(int[] core, int[] anotherCore) {
        for (int i = 0; i < core.length; i++) {
            try {
                //trying unisgned comparison to use all bits
                if (core[i] == anotherCore[i]) {
                    //keep going
                } else if ((core[i] < anotherCore[i]) ^ (core[i] < 0) ^ (anotherCore[i] < 0)) {
                    return 1;
                } else {
                    return -1;
                }
            } catch (IndexOutOfBoundsException e) {
                System.err.println(e.getMessage());
            }
        }
        return 0;
    }

    public static int computeHash(long[] core) {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(core);
        return hash;
    }

    public static CharSequence decodeCore(int encodedSequenceLength, int ecodedCore) {
        StringBuilder sb = new StringBuilder();
        int startPrintingBitsFrom = encodedSequenceLength * 2 - 1;
//        if (kmerCoreBits == 170160154469L) {
//            String toBinaryString = Long.toBinaryString(kmerCoreBits);
//            System.err.println(toBinaryString);
//            for (int i = 64; i > -1; i--) {
//                System.err.printf("%2d %d\n",i,(kmerCoreBits >>> i & 1));
//            }
//            System.err.println("DONE");
//        }

        for (int j = startPrintingBitsFrom; j > -1; j -= 2) {
            int twoBits = ecodedCore >> j - 1 & 3;
            switch (twoBits) {
                case 0:
                    sb.append("A");
                    break;
                case 1:
                    sb.append("C");
                    break;
                case 2:
                    sb.append("G");
                    break;
                case 3:
                    sb.append("T");
                    break;
                default:
                    break;
            }
        }
        return sb;
    }

    public static int encodeCore(CharSequence sequence) {
//        System.err.println("Encoding seq len=" + kmerString.length());
        int stringLength = sequence.length();
//        int currentInt = 0;
        int position = 0;
        int kmerCoreBits = 0;
        while (position < stringLength) {
            kmerCoreBits <<= 1;
            switch (sequence.charAt(position)) {
                case 'A':
                case 'a':
                    //if A : 00
                    kmerCoreBits <<= 1;
                    break;
                case 'C':
                case 'c':
                    //if C : 01
                    kmerCoreBits <<= 1;
                    kmerCoreBits++;
                    break;
                case 'G':
                case 'g':
                    //if G : 10
                    kmerCoreBits++;
                    kmerCoreBits <<= 1;
                    break;
                case 'T':
                case 't':
                    //if T : 11
                    kmerCoreBits++;
                    kmerCoreBits <<= 1;
                    kmerCoreBits++;
                    break;
                default:
                    System.err.println("Failed ecoding kmerstring to long....");
                    System.err.println("Offending char: " + sequence.charAt(position));
                    System.err.println("in " + sequence);
                    System.err.println("....exiting");
                    System.exit(1);
            }
            position++;
        }
//        String decodeCore = decodeCore(stringLength);
//        if (!decodeCore.equals(coreString)) {
//            System.err.println("Error encoding/decoding " + kmerCoreBits);
//            System.err.println(coreString + " <-core");
//            System.err.println(decodeCore + " <-decoded");
//            decodeCore = decodeCore(stringLength);
//        }
        return kmerCoreBits;
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet();
        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
//        optSet.addOpt(new Opt('A', "expected-adapter", "Expected adapter sequence (a fragment will do)", 1).setDefaultValue("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"));
//        optSet.addOpt(new Opt(null, "adapter-prefix-length", "Length of the adapter prefix used to identify 3' read-through", 1).setDefaultValue(9));
//        optSet.addOpt(new Opt('B', "blank-samples-name", "Name denoting blank samples in the key file. Name will by extended with remaining key-file fields", 1).setDefaultValue("Blank"));
//        optSet.addOpt(new Opt('U', "in-buffer-size", "Number of FASTQ records (reads or pairs depending on input) "
//            + "passed to in-queue", 1024, 128, 8092));
//        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up",
//            2, 1, 256));
//        //TRIMMING AND LENGTH
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Trimming and length settings]");
//        optSet.addOpt(new Opt('b', "keep-barcodes", "Do not trim barcodes"));
//        optSet.addOpt(new Opt('a', "keep-adapters", "Do not trim adapters found next to PstI and MspI sites "));
//        optSet.addOpt(new Opt('p', "keep-non-PstI-starting", "Keep reads (or pairs) which do not start with 'barcodeTGCAG'"));
//
//        int footId = 1;
//        String footText1 = "Note that certain combinations of min-length-* settings can lead to both mates of a pair ending up in SE/orphans output file.";
        optSet.addOpt(new Opt(null, "prefix-len", "Will determine the number of separate arrays to hold mers", 2, 1, 15, 1, 1));
        optSet.addOpt(new Opt(null, "postfix-len", "Will determine the number slots per array", 2, 1, 15, 1, 1));
//        optSet.addOpt(new Opt('e', "min-length-pair-each", "Only output a read pair if length of each is no less than <arg> bp, otherwise process as single", 1, 1, null, 1, 1).addFootnote(footId, footText1));
//        optSet.addOpt(new Opt('s', "min-length-pair-sum", "Only output a read pair if combined length is no less than <arg> bp, otherwise process as single", 2, 2, null, 1, 1).addFootnote(footId, footText1));
//        //RUNTIME
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime settings]");
//        optSet.addOpt(new Opt('c', "only-count", "Do not output reads"));
        optSet.addOpt(new Opt('t', "threads", "Number of threads", 1, 1, Runtime.getRuntime().availableProcessors(), 1, 1));
//        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
//
//        //OUTPUT
//        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Output settings]");
//        optSet.addOpt(new Opt('o', "out-dir", "Output directory", 1).setDefaultValue("out_split"));
//        optSet.addOpt(new Opt('x', "out-suffix-r1", "Output file suffix for R1 reads", 1).setDefaultValue("_R1.fastq.gz"));
//        optSet.addOpt(new Opt('X', "out-suffix-r2", "Output file suffix for R2 reads", 1).setDefaultValue("_R2.fastq.gz"));
//        optSet.addOpt(new Opt('S', "out-suffix-se", "Output file suffix for SE/orphaned reads", 1).setDefaultValue("_SE.fastq.gz"));
//        optSet.addOpt(new Opt(null, "append", "If output file(s) exist(s) for a given sample, append"));
//        optSet.addOpt(new Opt(null, "force", "If output file(s) exist(s) for a given sample, force overwrite"));
//        optSet.addOpt(new Opt('M', "[TODO] matchless-output", "Output reads with unmatched barcodes to R1/R2/SE file(s) prefixed with <arg>. If not set, these reads will be discarded", 1));
//        footId++;
//        String footText2 = "Consider increasing to sacrifice memory for speed. Decrease if encountering 'out of memory' errors.";
//        optSet.addOpt(new Opt('u', "out-buffer-size", "Number of FASTQ records (reads or pairs) "
//            + "passed to out-queue", 1024, 64, 8092).addFootnote(footId, footText2));
//        optSet.addOpt(new Opt('q', "out-queue-capacity", "Maximum number of buffers put on queue for writer threads to pick-up",
//            2, 1, 32).addFootnote(footId, footText2));
//
//        //POSITIONAL
        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAMEs", "names of input files", 1, (int) Short.MAX_VALUE));
        return optSet;
    }

//    private void addMer(PairMer[][] arr, int prefixLen, int postfixLen, String coreReminder, char clip, boolean left, Long[] stats) {
//        int prefix = encodeCore(coreReminder.subSequence(0, prefixLen));
//        int prefixp = encodeCore(coreReminder.subSequence(prefixLen, prefixLen + postfixLen));
//
//        find:
//        if (arr[prefix][prefixp] == null) {
//            try {
////                    arr[prefix][prefixp] = new PairMerIntArrEncoded(inputKmerSequence, prefixLen + postfixLen, inputKmerSequence.length()-1, false, 1);                    
//                arr[prefix][prefixp] = new NewPairMerIntArrEncoded(coreReminder, clip, left, 1);
//                stats[0]++;
//            } catch (NonACGTException ex) {
//                ex.printStackTrace();
//            }
//        } else {
//            PairMer previous = arr[prefix][prefixp];
//            while ((previous = previous.getNextPairMer()) != null) {
//                if (previous.decodeCore(coreReminder.length()).equals(coreReminder)) {
//                    //same core - TODO - add PM to core                    
//                    stats[1]++;
//                    break find;
//                } else {
////                        System.err.println(previous.decodeCore(coreReminder.length()) + " != " + coreReminder);
//                }
//            } 
//            try {
//                //replace first elem in linked list
////                PairMer stored = arr[prefix][prefixp];                
////                NewPairMerIntArrEncodedWithRef toStore = new NewPairMerIntArrEncodedWithRef(coreReminder, clip, left, 1, arr[prefix][prefixp]);
//                arr[prefix][prefixp] = new NewPairMerIntArrEncodedWithRef(coreReminder, clip, left, 1, arr[prefix][prefixp]);
//                stats[2]++;
//            } catch (NonACGTException ex) {
//                ex.printStackTrace();
//            }
//        }
//    }

    public CoreCoder(String[] args, String callerName, String toolName) {
        String TOOL_NAME = callerName + " " + toolName;
        int HELP_WIDTH = 150;
        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
//        readArgValues(optSet);

        int prefixLen = (int) optSet.getOpt("prefix-len").getValueOrDefault();;
        int prefixMax = (int) Math.pow(4, prefixLen);
        int postfixLen = (int) optSet.getOpt("postfix-len").getValueOrDefault();;
        int postfixpMax = (int) Math.pow(4, postfixLen);
        
        ArrayList<String> inputFileNamesList = new ArrayList<>();
        ArrayList<PositionalOpt> positionalOptsList = optSet.getPositionalOptsList();
        for (PositionalOpt po : positionalOptsList) {
            if (po.getValues() != null) {
                inputFileNamesList.addAll(po.getValues());
            }
        }

        NewPairMerArrays arr = new NewPairMerArrays(prefixLen, postfixLen, TOOL_NAME);
        Random random = new Random();
//        Long count = 0L;
//        Long dups = 0L;
//        Long diffs = 0L;
//        Long stats[] = new Long[]{0L,0L,0L};
            AtomicLongArray stats = new AtomicLongArray(3);




int INPUT_QUEUE_SIZE = 2;

 BlockingQueue inputQueue = new ArrayBlockingQueue(INPUT_QUEUE_SIZE);

        try {
            //READ INPUT AND POPULATE PairMers MAP
            int threads = (int)optSet.getOpt("threads").getValueOrDefault();
            int INPUT_BUFFER_SIZE = 1024;
            int MIN_KMER_FREQUENCY =1 ;
            Reporter.report("[INFO]", "Allocated " + threads + " thread(s) to map populating", TOOL_NAME);
            ArrayList<Future<?>> futures = new ArrayList<>(threads + 1);
            final ExecutorService readAndPopulateExecutor = new ThreadPoolExecutor(threads + 1, threads + 1, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

            //SPAWN INPUT READING THREAD
            ArrayList<Integer> kSizes = new ArrayList<>();
            InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, kSizes, inputFileNamesList, INPUT_BUFFER_SIZE, TOOL_NAME);

            Future<?> future = readAndPopulateExecutor.submit(inputReaderProducer);
            futures.add(future);
            boolean splitInputSequenceIntoKmers = true;

            //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
            long timeStart = System.currentTimeMillis();
            int count = 0;
            while (inputReaderProducer.getGuessedInputFormat() == null) {
                try {
                    //IF nothing happens after 5 seconds
                    if (System.currentTimeMillis() - timeStart > 5000 && (count++ % 50 == 0)) {
                        Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", TOOL_NAME);
                    }
                    Thread.sleep(100); //wait for 1/10 of a second
                } catch (InterruptedException ex) {
                }
            }
//            Integer kSizeFromInput = null;
            if (inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.InFormat.KMERS)) {
                splitInputSequenceIntoKmers = false;
            }

            //SPAWN THREADS TO POPULATE MAP
            for (int i = 0; i < threads; i++) {
//                PairMerMapPopulatorConsumer consumer = new PairMerMapPopulatorConsumer(inputQueue, pairMerMaps,
//                        splitInputSequenceIntoKmers, kSizes, MIN_KMER_FREQUENCY, false);
                NewPairMerMapPopulatorConsumer consumer = new NewPairMerMapPopulatorConsumer(inputQueue, arr, false, kSizes, i, false, stats, prefixLen, postfixLen);
                futures.add(readAndPopulateExecutor.submit(consumer));
            }















            readAndPopulateExecutor.shutdown();
            try {
                for (Future<?> f : futures) {
                    f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                }
                readAndPopulateExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Reporter.report("[ERROR]", "PairMerSet populator interrupted exception!", TOOL_NAME);
            } catch (ExecutionException ex) {
                Reporter.report("[ERROR]", "PairMerSet populator execution exception!", TOOL_NAME);
                ex.printStackTrace();

            } catch (TimeoutException ex) {
                Logger.getLogger(KmerExtender.class.getName()).log(Level.SEVERE, null, ex);
                Reporter.report("[ERROR]", "PairMerSet populator timeout exception!", TOOL_NAME);
            }

//            if (pairMersMap.isOutOfMemory()) {
//                Reporter.report("[ERROR]", "Terminating, out of memory while populating the Map, k=" + KMER_LENGTH + ", n=" + NumberFormat.getIntegerInstance().format(pairMersMap.getPairMersSkipListMap().size()), TOOL_NAME);
//                System.exit(1);
//            }
        } catch (OutOfMemoryError e) {
            Reporter.report("[ERROR]", "Out of memory error!", TOOL_NAME);
//            printStatus("Reference set n = " + NumberFormat.getNumberInstance().format(clipMersMap.getClipMersMap().size()));
            System.exit(1);
        }








        
//        int randomInts = prefixMax * postfixpMax*100;
//        for (int i = 0; i < randomInts; i++) {
////            CharSequence kmerSequence = decodeCore(16, random.nextInt(1024));
//            CharSequence inputKmerSequence = decodeCore(16, random.nextInt());
//
//            arr.addKmer(inputKmerSequence, stats);
//
//            //elem not seen before, store
//        }
//        Reporter.report("[INFO]", randomInts + " input kmers", TOOL_NAME);
        Reporter.report("[INFO]", NumberFormat.getIntegerInstance().format(stats.get(0)) + " PairMers loaded without conflict" , TOOL_NAME);
        Reporter.report("[INFO]", NumberFormat.getIntegerInstance().format(stats.get(1)) + " duplicated PairMer cores", TOOL_NAME);
        Reporter.report("[INFO]", NumberFormat.getIntegerInstance().format(stats.get(2)) + " distinct PairMer cores at given prefix", TOOL_NAME);
        Reporter.report("[INFO]", NumberFormat.getIntegerInstance().format((stats.get(0)) + stats.get(1) + stats.get(2))+" total PairMers ", TOOL_NAME);
//        Reporter.report("[INFO]", "Load level 1 = " +  perc.format(((double) lev1load) / arr.length), TOOL_NAME);

        arr.printStats();
    }

//    private class PMO {
//
//        private final int suffix;
//
//        public PMO(int suffix) {
//            this.suffix = suffix;
//        }
//
//        public int getSuffix() {
//            return suffix;
//        }
//
//        public PMO getNextPairMer() {
//            return null;
//        }
//
//    }
//
//    private class PMO1 extends PMO {
//
//        private final PMO next;
//
//        public PMO1(int suffix, PMO next) {
//            super(suffix);
//            this.next = next;
//        }
//
//        @Override
//        public PMO getNextPairMer() {
//            return next;
//        }
//
//    }
}
