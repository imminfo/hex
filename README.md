[![Licence](https://img.shields.io/hexpm/l/plug.svg)](http://www.apache.org/licenses/LICENSE-2.0)

# HEX
HEX(HLA Extractor) is a simple HLA typing tool. It is designed to ???


## Installation and dependencies
It depends on 
- Python3
- Perl
- bwa-mem
- ncbi-blast


## Workflow

### Step 1. Building the database for HLA
`hex_build.py`

By default with the HEX comes the preprocessed database for IMGT <name here.gz> from <date here>.

#### 1a. Loading the IMGT database with HLA nucleotide sequences

The script loads IMGR database files from GitHub bu default. You can specify URLs with alleles by yourself (see the `hlalinks.txt` file).

1. HLA 1st type -> 2-3-4 exons only plus a little bit of the 1st exon
2. HLA 2nd type -> 1-2-3 exons only

#### 1b. Construction of the C-value matrix for alleles

C-value ("confuse value") is the probability that one allele is replaced with another by an error.

estimating 1st and 5th exons for missing data


### Step 2. HLA extraction
`hex_extract.py`

#### 2a. Preprocessing input .fastq files

1. Align sequences using `bwa-mem` to reference genes.
2. Split initial .fastq files to files related to a specific genes.

#### 2b. Make consensuses

1. Make consensuses sequences from each gene-specific file.

#### 2c. Generate HLA candidates

1. Align candidate sequences to the HLA database using `BLAST`.
