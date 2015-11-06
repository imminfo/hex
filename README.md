[![Licence](https://img.shields.io/hexpm/l/plug.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0)

# HEX
HEX (HLA Extractor) is a HLA typing software tool designed for HLA typing NGS data obtained with the novel developed technology. HEX was developed in the [Laboratory of Comparative and Functional Genomics](http://labcfg.ibch.ru/lcfg.html).

## Installation and dependencies
It depends on:
- [Python 3](https://www.python.org/downloads/) is the the main programming language of the tool.
- [bwa-mem](http://bio-bwa.sourceforge.net) is used for mapping the input sequences to the reference HLA sequences.
- [Perl](https://www.perl.org/get.html) is used in processing the `bwa-mem` mapping results.
- [ncbi-blast](https://www.ncbi.nlm.nih.gov/books/NBK279671/) is used for aligning the resulting HLA candidates to the alleles.

## Run HEX

HEX already comes with prepocessed database with alleles and reference sequences. To start HEX you just need to run:

```
$./hex_extract <path to a folder with raw .fastq files>
```

## Algorithm

### Step 1. Allele database preprocessing

In order to build the database with reference sequences HEX need to download the IMGT HLA data and process it. To do so you just need to run:

```
$./hex_build [-i <optional input file with URLs>] [-o <optional output folder>]
```

`-i` - a path to a file with links to files with HLA (default - `hlalinks.txt`).
`-o` - output folder, to which the database will be written.

### Step 2. HLA extraction

#### 2a. Preprocess input FASTQ files

For each input sample do:

1. Map sequences using `bwa-mem` to the reference sequences.

2. For each sequence fill the gaps with `Xs` and write to the file related to the specific HLA class.

#### 2b. Consensus sequences construction

Make consensuses sequences from each HLA class-specific file:

1. Choose the non-conservative position with the best coverage.

2. For each direction (i.e., two) make the two best sequences using the markov chains logic.

3. Merge the corresponding consensus sequences by the start nucleotide.

#### 2c. Generate HLA candidates

1. Align consensus sequences to the HLA database using `BLAST`.

## Input NGS data

The input NGS data is foollows:

  ???
  
  ???
  
  ???

