- [`yakat` - yet another k-mer analysis toolkit?](#yakat---yet-another-k-mer-analysis-toolkit)
- [Dependencies](#dependencies)
- [You may also need](#you-may-also-need)
- [Compile and build](#compile-and-build)
- [Run](#run)
- [Usage in pipelines](#usage-in-pipelines)
- [Development](#development)
- [Bugs, cryptic errors, general enquires](#bugs-cryptic-errors-general-enquires)
# `yakat` - yet another  k-mer analysis toolkit?

Or is it?
Perhaps better described as _not quite_ or _not just_ another k-mer analysis toolkit.
This software included several k-mer related tools, mostly trying to complement or functionally extend existing tools.
Written in Java these are not necessarily most memory frugal, but most are pretty fast, thanks to leveraging producer-consumer approach to multi-threading in Java.
So if you happen to have access to a fat node with lots of CPUs as I did at some point, you may find these quite handy.
In addition to k-mer based tools, there are several other modules providing efficient, multi-threaded solutions for common bioinformatics task such as read id matching/filtering.

# Dependencies

You'll need Java 8 with [ant](https://ant.apache.org/) for compiling and building, or just Java for running the pre-compiled binary.
Normal usage is in Linux environment but most modules should work on other systems with no or few adjustments required, such as explicitly specifying input and output files other than `/dev/stdin` and `/dev/stdout`.

# You may also need

A k-mer counter. Go for [KMC](https://github.com/refresh-bio/KMC).
Some yakat modules will k-merize input as needed, but use KMC or another, dedicated k-mer counter whenever a set of k-mers is a desired input, especially for the `kextend` module.

# Compile and build

After cloning or downloading this repository, run `ant compile && ant jar`, this should generate the self-contained Java executable `dist/yakat.jar`.

# Run

```
java -jar dist/yakat.jar

Usage: java -jar yakat.jar <module>
k-mer based modules
  kextend       : extend k-mers to unambiguous contigs or extend input "seed" sequences only
  snpmers       : given parental SNPs (e.g. from LNISKS), corresponding FASTA sequences and sets of k-mers,
                  call offspring genotypes by overlapping their k-mers with parental SNP sequences
  freqmers      : given a set of sequences and set(s) of k-mers
                  report k-mer coverage and frequency for the input sequences
  kmatch        : match/filter/bait FAST(A|Q) sequences based on contained k-mers (or lack thereof)
  seedmers      : [PROTOTYPE] given seed seequences interrogare sets of k-mers
                  to genotype presumed mutations at positions k bases from the seed edges

FASTQ processing modules:
  idmatch       : match FASTQ records by id
  split         : split FASTQ GBS/ddRAD reads by barcodes, trim barcodes and adapters

MPILEUP processing modules:
  pileupstats   : extract some stats from (m)pileup
  pmpileup      : count and call bases from (m)pileup

HMMER, VSEARCH output post-processing modules:
  hmmerdoms     : Process and group domain-hits from HMMER
  vclusters     : call variants from VSEARCH clustering msa output

Miscellaneous modules:
  allorfs       : Identify and extract all (longest) ORFs from a genome
  kexpress      : [UNTESTED] Calculate expresion values (TPM) from read counts
                    and (not quite) GFF description of features

Info:
  version       : print the version and exit
```

As you can see, not everything k-mer based, additional, deprecated modules are still lurking in the code.

The first time around you may want to run a specific module do it with `-h` or `--help` as some modules kick off by reading from `/dev/stdin`.


```
java -jar dist/yakat.jar kmatch -h
```

# Usage in pipelines

Note that due to specific modules' original application within larger pipelines they occasionally expect slightly modified file formats.
For example by FASTQ record we mean a single read (pair of reads) whose four (eight) lines have been placed on a single line using tab as a separator.

To filter (in/out) reads based on matching k-mers you may run the `kmatch` module like this:

```
zcat reads.fastq.gz \
  | paste - - - - \
  | java -jar scripts/yakat.jar kmatch \
    --k-mers 4_reads.fastq \
    --k-mer-length 50 \
  | tr '\t' '\n' > filtered.fastq
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     Start populating k-mers set [1.92 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     Input format guessed: FASTQ [2.4 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     4 FASTQ read-in [2.4 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     Finished populating k-mers set, n=204 [2.88 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     Input format guessed: FASTQ_SE_ONE_LINE [3.36 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     2,048 FASTQ_SE_ONE_LINE read-in so far [4.8 MB]
2018-12-12 Wed 14:21:27 [yakat kmatch]         [INFO]     4,096 FASTQ_SE_ONE_LINE read-in so far [7.21 MB]
2018-12-12 Wed 14:21:28 [yakat kmatch]         [INFO]     8,192 FASTQ_SE_ONE_LINE read-in so far [21.74 MB]
2018-12-12 Wed 14:21:29 [yakat kmatch]         [INFO]     12,292 FASTQ_SE_ONE_LINE read-in [55.32 MB]
2018-12-12 Wed 14:21:30 [yakat kmatch]         [INFO]     Finished outputing records, n=17 [92.05 MB]
```

The content of `4_reads.fastq` was k-merized and used for matching reads streamed from stdin.

To filter in/out pairs of reads you could:

```
paste <(zcat R1.fq.gz | paste - - - - ) \
      <(zcat R2.fq.gz | paste - - - - ) \
  | java -jar scripts/yakat.jar kmatch \
    --k-mers 4_reads.fastq \
    --k-mer-length 50 \
  | tee >(cut -f 1-4 -d$'\t' | tr '\t' '\n' > filtered_R1.fq) \
  | cut -f 5-8 -d$'\t' | tr '\t' '\n' > filtered_R2.fq \
```

# Development

Developed as a NetBeans project and can be loaded as such. Stress-tested on ubuntu 14.04 with Java 1.8 build 74.

# Bugs, cryptic errors, general enquires

[Submit an issue](https://github.com/rsuchecki/yakat/issues/new)




