[![DOI](https://zenodo.org/badge/161127471.svg)](https://zenodo.org/badge/latestdoi/161127471)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/rsuchecki/yakat/CI)](https://github.com/rsuchecki/yakat/actions)
[![Latest GitHub tag](https://img.shields.io/github/tag/rsuchecki/yakat.svg?label=latest%20release&logo=github)](https://github.com/rsuchecki/yakat/releases)
[![GitHub commits since latest tag/release](https://img.shields.io/github/commits-since/rsuchecki/yakat/latest.svg?logo=github)](https://github.com/rsuchecki/yakat/releases)


# Table of Contents <!-- omit in toc -->

- [`yakat` - yet another  k-mer analysis toolkit?](#yakat---yet-another--k-mer-analysis-toolkit)
- [Getting started](#getting-started)
  - [Dependencies](#dependencies)
  - [Compile and build and run](#compile-and-build-and-run)
- [Get cracking](#get-cracking)
  - [You may also need](#you-may-also-need)
  - [Selected use cases](#selected-use-cases)
    - [`freqmers`](#freqmers)
    - [`kextender` - default mode](#kextender---default-mode)
    - [`kextender` - FASTA seed extension mode](#kextender---fasta-seed-extension-mode)
    - [`kmatcher` - see usage in pipelines](#kmatcher---see-usage-in-pipelines)
    - [`snpmers`](#snpmers)
    - [`vclusters`](#vclusters)
  - [Usage in pipelines](#usage-in-pipelines)
- [Development](#development)
- [Bugs, cryptic errors, general enquires, existential angst?](#bugs-cryptic-errors-general-enquires-existential-angst)
# `yakat` - yet another  k-mer analysis toolkit?

Or is it?
Perhaps better described as _not quite_ or _not just_ another k-mer analysis toolkit.
This software includes several k-mer related tools, mostly trying to complement or functionally extend existing tools.
Written in Java these are not necessarily most memory frugal, but most are fast, thanks to leveraging producer-consumer approach to multi-threading in Java.
So if you happen to have access to a fat node with lots of CPUs as I did at some point, you may find these quite handy.
In addition to k-mer based tools, there are several other modules providing efficient, multi-threaded solutions for common bioinformatics task such as read id matching/filtering.

# Getting started

## Dependencies

You'll need Java 8 with [ant](https://ant.apache.org/) for compiling and building, or just Java for running the pre-compiled binary.
Normal usage is in Linux environment but most modules should work on other systems with no or few adjustments required, such as explicitly specifying input and output files other than `/dev/stdin` and `/dev/stdout`.


## Compile and build and run

After cloning or downloading this repository, run

`ant jar` 

or, if that fails,  `ant -Dplatforms.JDK_8.home=${JAVA_HOME} jar` 

This should generate the Java executable `dist/yakat.jar`
and a self contained linux executable `yakat`

You can now either run `./yakat` or `java -jar dist/yakat.jar`.

In the former case you can pass Java VM options using `--JVM "<options>"`, for example: `./yakat --JVM "-Xmx2G -Xms100m"`

You should see the following summary of available modules.

```
Usage:
  either:   yakat <module>
  or:       java -jar yakat.jar <module>

k-mer based modules
  freqmers      : given a set of sequences and set(s) of k-mers
                  report k-mer coverage and frequency for the input sequences
  kextend       : extend k-mers to unambiguous contigs or extend input "seed" sequences only
  kmatch        : match/filter/bait FAST(A|Q) sequences based on contained k-mers (or lack thereof)
  seedmers      : [PROTOTYPE] given seed seequences interrogare sets of k-mers
                  to genotype presumed mutations at positions k bases from the seed edges
  snpmers       : given parental SNPs (e.g. from LNISKS), corresponding FASTA sequences and sets of k-mers,
                  call offspring genotypes by overlapping their k-mers with parental SNP sequences

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

As you can see, not everything is k-mer based, additional, deprecated modules are still lurking in the code as well.

The first time you run a specific module do it with `-h` or `--help` as some modules kick off by reading from `/dev/stdin`.

# Get cracking

## You may also need

A k-mer counter. Go for [KMC](https://github.com/refresh-bio/KMC).
Some yakat modules will k-merize input as needed, but use KMC or another, dedicated k-mer counter whenever a set of k-mers is the desired input, especially for the `kextend` module.

## Selected use cases

### `freqmers`

This module can be used to quantify how sets of k-mers relate to / overlap with 
the set of input (FASTA) sequences. 

* Take FASTA sequences 
* Take a set of k-mers per sample of interest. 
* Record frequencies of k-mers overlapping the sequences 
* Report k-mer coverage and frequencies along the sequences 

### `kextender` - default mode

Among the available modules `kextender` is by far the most mature, if you are after no-nonsense, fast generation of unitigs from a set of Illumina reads all you need to do is:

* k-merize your reads with KMC (other k-mer counters are available)
* determine k-mer frequency cutoff `${MIN_FREQ}` to exclude low frequency k-mers which are likely error-induced by looking at the output of `kmc_tools histogram`
* pipe your k-mers from KMC database to `yakat kextend`

```sh
kmc_dump -ci${MIN_FREQ} db_basename /dev/stdout \
  | ./yakat kextend > unitigs
```

Note that the default output is one unitig (sequence) per line,
for FASTA output use `--fasta-out` flag.
You may also use `--min-length` to set the minimum length (bp) of an output unitig.


### `kextender` - FASTA seed extension mode

The `kextender` module can also extend existing _seed_ sequences. Note that such extensions never lead to the seed sequences being joined. 
The idea is to be able to unambiguously extend existing sequences - think extending with a unitig. 
Unlike the default mode which uses a single _k_ value in the extension process, the seed extension can make use of multiple 
values of _k_. Extensions of each end of a seed are computed for a range of values of _k_ and the longest of those is "attached" 
to the seed end. Exploring a large range of k values for a significant input [k-mers/FAST[A|Q]] will make the memory requirements explode, 
as it is done in parallel 

* in the interest of speed 
* to avoid excessive I/O 
* to preserve stdin handling

Example usage:

```sh
kmc_dump -ci${MIN_FREQ} db_basename /dev/stdout \
  | ./yakat kextend \
    --seed-file existing_sequences.fasta \
    --k-mer-min 20 \
    --k-mer-max 80 \
    --k-mer-step 1 \
    > extended_seeds.fasta 
```

### `kmatcher` - see [usage in pipelines](#usage-in-pipelines)

### `snpmers`

See [this example](https://github.com/rsuchecki/LNISKS/tree/d165629ae1d40ae7158819149d9b56fbba42e995#call-prioritization).


### `vclusters`

This module parses the MSA output of VSEARCH clustering
and calls variants within each cluster.



## Usage in pipelines

Note that due to specific modules' original application within larger pipelines they occasionally expect/produce slightly modified versions of common file formats.
Fear not, these are not whimsical modifications of accepted standards but rather alternative presentation of existing formats which facilitates parallelised processing and use of linux pipes rather than intermediary files. For example,

* by a FASTQ_SE_ONE_LINE record we mean a single read whose four lines have been placed on a single line using tab as a separator.
* by a FASTQ_PE_ONE_LINE record we mean a pair of reads whose eight lines have been placed on a single line using tab as a separator.

The idea is to wrap/unwrap these on the fly, either on the command line or using wrapper scripts.

To filter (in/out) reads based on matching k-mers you may run the `kmatch` module like this:

```
zcat reads.fastq.gz \
  | paste - - - - \
  | ./yakat kmatch \
    --k-mers 4_reads.fastq \
    --k-mer-length 50 \
  | tr '\t' '\n' > filtered.fastq
```

* Each read is unwrapped into the one-line-per-record form, processed using `kmatch` and surviving reads are unwrapped by replacing tabs with newline characters.
* The content of `4_reads.fastq` is k-merized and used for matching reads streamed from stdin.
* In practice the output could be piped into another process or perhaps a parallelized compression tool such as `pigz`.

Analogous operation for paired-end data:

```
paste <(zcat R1.fq.gz | paste - - - - ) \
      <(zcat R2.fq.gz | paste - - - - ) \
  | ./yakat kmatch \
    --k-mers 4_reads.fastq \
    --k-mer-length 50 \
  | tee >(cut -f 1-4 -d$'\t' | tr '\t' '\n' > filtered_R1.fq) \
  | cut -f 5-8 -d$'\t' | tr '\t' '\n' > filtered_R2.fq
```

# Development

Developed as a NetBeans project and can be loaded as such. Stress-tested on ubuntu 14.04 with Java 1.8 build 74.

# Bugs, cryptic errors, general enquires, existential angst?

[Submit an issue](https://github.com/rsuchecki/yakat/issues/new)




