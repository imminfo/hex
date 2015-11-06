[![Licence](https://img.shields.io/hexpm/l/plug.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0)

# HEX
HEX (HLA Extractor) is a HLA typing software tool designed for HLA typing NGS data obtained with the novel developed technology.

## Installation and dependencies
It depends on:
- Python3
- Perl
- bwa-mem
- ncbi-blast

## Run HEX

HEX already comes with prepocessed database with alleles and reference sequences. To start HEX you just need to run:

```
$./hex_extract <path to a folder with raw .fastq files>
```

## Algorithm

### Step 1. Allele database preprocessing
`$./hex_build [-i <optional input file with URLs>] [-o <optional output folder>]`

By default with the HEX comes the preprocessed database for IMGT <name here.gz> from <date here>.

#### 1a. Loading the IMGT database with HLA nucleotide sequences

The script loads IMGR database files from GitHub bu default. You can specify URLs with alleles by yourself (see the `hlalinks.txt` file).

1. HLA 1st type -> 2-3-4 exons only plus a little bit of the 1st exon
2. HLA 2nd type -> 1-2-3 exons only

### Step 2. HLA extraction

#### 2a. Preprocessing input .fastq files

1. Align sequences using `bwa-mem` to reference genes.
2. Split initial .fastq files to files related to a specific genes.

#### 2b. Make consensuses

1. Make consensuses sequences from each gene-specific file.

#### 2c. Generate HLA candidates

1. Align candidate sequences to the HLA database using `BLAST`.

## Input NGS data
