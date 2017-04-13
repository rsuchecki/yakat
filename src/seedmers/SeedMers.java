/*
 * Copyright 2017 Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>.
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
package seedmers;

import argparser.ArgParser;
import argparser.Opt;
import argparser.OptSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.zip.GZIPInputStream;
import shared.InputReaderProducer;
import shared.Message;
import shared.Reporter;
import shared.SequenceOps;
import shared.StdRedirect;

/**
 * Given a seed sequence, generate possible alternative sequences (altseeds)
 * such that a base x at position k bases from a sequence end is replaced with
 * base y \in {A,C,G,T}, y != x;
 *
 * Thread the input set(s) of k-mers through altseeds to genotype the respective source
 *
 * For now: focus on seeds length=2k-1, only focus on SNPs, ignore indels
 *
 * @author Radoslaw Suchecki <radoslaw.suchecki@adelaide.edu.au>
 */
public class SeedMers {

    private final String TOOL_NAME;
    private final int HELP_WIDTH = 200;
    private final int READER_BUFFER_SIZE = 8192;
    private final int WRITER_BUFFER_SIZE = 8192;
    private final int REPORTING_SHIFT = 6; //input progress reporting every READER_BUFFER_SIZE*REPORTING_FACTOR records
    private HashMap<String, ArrayList<AltSeedLink>> map;
    private ArrayList<Seed> seeds;

    public SeedMers(String[] args, String callerName, String toolName) {
        TOOL_NAME = callerName + " " + toolName;

        OptSet optSet = populateOptSet();
        ArgParser argParser = new ArgParser();
        argParser.processArgs(args, optSet, true, callerName, HELP_WIDTH);
        new StdRedirect(optSet, TOOL_NAME);
        if (optSet.getOpt("P").isUsed()) {
            optSet.printUserSettings(TOOL_NAME);
        }

        buildSeedMersMap(optSet);
//        removeSnpMersPoorlyCoveredInParents(optSet); //DONE AFTER FALSE POSITIVE FILTERING

        int minTotal = (int) optSet.getOpt("min-k-mer-frequency-sum").getValueOrDefault();
        int minMinor = (int) optSet.getOpt("min-k-mer-frequency-minor").getValueOrDefault();
        double minCoverage = (double) optSet.getOpt("min-snp-coverage").getValueOrDefault();
        double maxError = (double) optSet.getOpt("max-coverage-error").getValueOrDefault();

//        ArrayList<String> sampleNames = threadKmersThroughMap(optSet, kmersFileNames, minTotal, minMinor, minCoverage, maxError);
        threadKmersThroughMap(optSet);
        Reporter.report("[INFO]", "Done!", TOOL_NAME);
    }

    private OptSet populateOptSet() {
        OptSet optSet = new OptSet("description here");

        //INPUT
        optSet.setListingGroupLabel("[Input settings]");
        optSet.addOpt(new Opt('k', "k-mer-length", "It must match the size of the k-mers used to query the niks-snps", 1).setRequired(true).setMinValue(3).setMaxValue(255));
        optSet.addOpt(new Opt('f', "fasta-seeds", "The FASTA file containing seeds which will be interrogated for variation", 1).setRequired(true));
        optSet.addOpt(new Opt('K', "per-sample-k-mers", "A set (or sets) of k-mers to be threaded through the map of k-mer-links to presumed SNPs").setMinValueArgs(1).setMaxValueArgs(Integer.MAX_VALUE));
        optSet.addOpt(new Opt('U', "in-buffer-size", "Size of buffers put on in-queue ", 8192, 128, 65535));
        optSet.addOpt(new Opt('Q', "in-queue-capacity", "Maximum number of buffers put on queue for processing threads to pick-up", 64, 1, 256));

        int footId = 1;
        String footText = "k-mer frequency corresponding to a SNP allele is obtained by taking a median of frequencies of all k-mers overlapping the underlying SNP";

        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Variant calling and reporting]");
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-sum", "Minimum frequency of k-mers which overlap "
                + "with a SNP (minor+major allele)", 1).setMinValue(1).setDefaultValue(5).addFootnote(footId, footText));
        optSet.addOpt(new Opt(null, "min-k-mer-frequency-minor", "Minimum frequency of k-mers "
                + "which support the minor allele", 1).setMinValue(1).setDefaultValue(2).addFootnote(footId, footText));
        optSet.addOpt(new Opt(null, "min-snp-coverage", "At least <arg> fraction of unique snp-covering k-mers must be present in a genotyped dataset (for each allele)", 1).setDefaultValue(0.9));
        optSet.addOpt(new Opt(null, "max-coverage-error", "Maximum fraction of unique snp-covering k-mers that will be ignored as errorneous assignment and not considered for genotyping", 1).setMinValue(0.00).setDefaultValue(0.05));
//
        optSet.setListingGroupLabel(optSet.incrementLisitngGroup(), "[Runtime and output settings]");
        optSet.addOpt(new Opt('o', "stdout-redirect", "Redirect stdout to this file", 1));
        optSet.addOpt(new Opt('e', "stderr-redirect", "Redirect stderr to this file", 1));
        optSet.addOpt(new Opt(null, "out-fasta", "Output relevant sequences to this file (in FASTA format)", 1));
        optSet.addOpt(new Opt(null, "out-calls", "Output calls to this file ", 1).setDefaultValue("/dev/stdout"));
        optSet.addOpt(new Opt(null, "out-calls-AB", "Output calls to this file (in AB format)", 1));
        optSet.addOpt(new Opt(null, "out-calls-IUPAC", "Output calls to this file (in IUPAC format) - due to the limitations "
                + "of this format indels will be lost", 1));
        optSet.addOpt(new Opt('P', "print-user-settings", "Print the list of user-settings to stderr and continue executing"));
        optSet.addOpt(new Opt('D', "debug", "Print additional info for debugging purposes"));
//        boolean positionalArgumentRequired = true;
//        optSet.addPositionalOpt(new PositionalOpt("INPUT_FILENAME", "name of input file ", 1, positionalArgumentRequired));
        return optSet;
    }

    /**
     * Read FASTA seeds, generate altseeds, generate k-mers from altseeds and
     * link those k-mers back to those altseeds (or just to the seed)
     *
     * @param optSet
     */
    private void buildSeedMersMap(OptSet optSet) {
        seeds = new ArrayList<>();
        map = new HashMap<>();
        String seedsFileName = (String) optSet.getOpt("f").getValueOrDefault();
        int k = (int) optSet.getOpt("k").getValueOrDefault();
        BufferedReader bufferdReader = null;
        try {
            String inputLine;
            if (seedsFileName.endsWith(".gz")) {
                InputStream gzipStream = new GZIPInputStream(new FileInputStream(seedsFileName), READER_BUFFER_SIZE);
                bufferdReader = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"), READER_BUFFER_SIZE);
            } else {
                bufferdReader = new BufferedReader(new FileReader(new File(seedsFileName)), READER_BUFFER_SIZE);
            }
//            Pattern spliPattern = Pattern.compile("\t");
            String id = null;
            StringBuilder seqBuilder = null;
            while ((inputLine = bufferdReader.readLine()) != null) {
                if (inputLine.startsWith(">")) {
                    if (id != null) {
                        Seed seed = new Seed(id, seqBuilder, k);
                        seeds.add(seed);
                        kmerizeAndLink(k, seed);
                    }
                    id = inputLine.trim();
                    seqBuilder = new StringBuilder();
                } else {
                    seqBuilder.append(inputLine);
                }
            }
            Seed seed = new Seed(id, seqBuilder, k);
            seeds.add(seed);
            kmerizeAndLink(k, seed);            
        } catch (FileNotFoundException ex) {
            Reporter.report("[ERROR]", "File not found: " + seedsFileName, TOOL_NAME);
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (bufferdReader != null) {
                    bufferdReader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        Reporter.report("[INFO]", "Seed-mers map populated, n="+NumberFormat.getIntegerInstance().format(map.size()), TOOL_NAME);
//        for (Map.Entry<String, ArrayList<AltSeedLink>> entry : map.entrySet()) {
//            String key = entry.getKey();
//            ArrayList<AltSeedLink> value = entry.getValue();
//            System.err.println(key);
//        }
    }

    /**
     * Consider multi-threading
     *
     * @param k
     * @param seed
     */
    private void kmerizeAndLink(int k, Seed seed) {
        CharSequence sequence = seed.getSequence();
        int snpSite = k - 1;
        int maxKmer = Math.min(sequence.length() - k + 1, snpSite + 1);
        int startAt = Math.max(0, snpSite - k + 1);

        Character altChars[] = new Character[]{'A', 'C', 'G', 'T'};
        for (Character base : altChars) {
            StringBuilder altSeq = new StringBuilder(sequence.subSequence(0, snpSite));
            altSeq.append(base);
            altSeq.append(sequence.subSequence(snpSite+1, sequence.length()));
            //k-merize
            for (int i = startAt; i < maxKmer; i++) {
                CharSequence kmer = altSeq.subSequence(i, i + k);
                String canonical = SequenceOps.getCanonical(kmer.toString());
                ArrayList<AltSeedLink> altSeeds = map.get(canonical);
                if (altSeeds == null) {
                    altSeeds = new ArrayList<>();
                }
                altSeeds.add(new AltSeedLink(seed, i, base, canonical.compareTo(kmer.toString()) == 0));
                map.put(canonical, altSeeds);
            }
        }
    }

    private ArrayList<String> threadKmersThroughMap(OptSet optSet) {
        ArrayList<String> samples = new ArrayList<>();
        int IN_Q_CAPACITY = (int) optSet.getOpt("Q").getValueOrDefault();
        int IN_BUFFER_SIZE = (int) optSet.getOpt("U").getValueOrDefault();
        ArrayList<String> kmersFileNames = (ArrayList<String>) optSet.getOpt("K").getValues();

        BlockingQueue inputQueue = new ArrayBlockingQueue(IN_Q_CAPACITY);
//            boolean stranded = false;
        int ioThreads = 1 + 1;
        ArrayList<Future<?>> ioFutures = new ArrayList<>(ioThreads);
        final ExecutorService ioExecutorService = new ThreadPoolExecutor(ioThreads, ioThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());

        //SPAWN INPUT READING THREAD
        InputReaderProducer inputReaderProducer = new InputReaderProducer(inputQueue, kmersFileNames, true, REPORTING_SHIFT, TOOL_NAME, IN_BUFFER_SIZE);
        ioFutures.add(ioExecutorService.submit(inputReaderProducer));

        long timeStart = System.currentTimeMillis();
        int count = 0;
        while (inputReaderProducer.getGuessedInputFormat() == null) {
            try {
                //IF nothing happens after 5 seconds
                if (System.currentTimeMillis() - timeStart > 2500 && (count++ % 25 == 0)) {
                    Reporter.report("[WARNING]", "Stuck or waiting for standard input stream...", getClass().getSimpleName());
                }
                Thread.sleep(100); //wait for 1/10 of a second
            } catch (InterruptedException ex) {
            }
        }
//        //ENSURING WE KNOW THE INPUT FORMAT BEFORE CONSUMER THREADS ARE SPAWNED
//        if (!inputReaderProducer.getGuessedInputFormat().equals(InputReaderProducer.InFormat.KMERS)) {
//            Reporter.report("[FATAL]", "Only k-mer sets accepted as input. Guessed format: " + inputReaderProducer.getGuessedInputFormat(), getClass().getSimpleName());
//        } else {
        //Start KmergerConsumerProducer and OutputWriterConsumer threads
        int threads = 1;
        final ExecutorService execService = new ThreadPoolExecutor(threads, threads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        ArrayList<Future<?>> futures = new ArrayList<>(threads);
//        AtomicInteger splitterThreads = new AtomicInteger(SPLITTER_THREADS);
        ArrayList<Message> finalMessages = new ArrayList<>(threads * 5);
//        for (int i = 0; i < threads; i++) {
        futures.add(execService.submit(new CallerConsumer(inputQueue, TOOL_NAME, samples, map, seeds)));
//        futures.add(execService.submit(new CallerConsumer(inputQueue, TOOL_NAME, samples, map, minTotal, minMinor, minCoverage, minError)));
//        }
        execService.shutdown();
        ioExecutorService.shutdown();

        try {

            for (Future<?> f : futures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            execService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            for (Future<?> f : ioFutures) {
                f.get(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }
            ioExecutorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Reporter.report("[ERROR]", "interrupted exception!", getClass().getSimpleName());
        } catch (ExecutionException ex) {
            Reporter.report("[ERROR]", "execution exception! " + ex.getCause().getMessage(), getClass().getSimpleName());
            ex.printStackTrace();
        } catch (TimeoutException ex) {
            Reporter.report("[ERROR]", "timeout exception!", getClass().getSimpleName());
        }
        for (Message fm : finalMessages) {
            Reporter.report(fm.getLevel().toString(), fm.getBody(), fm.getCaller());
        }

        return samples;
    }

}
