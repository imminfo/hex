[![Licence](https://img.shields.io/hexpm/l/plug.svg)](http://www.apache.org/licenses/LICENSE-2.0)

# HEX
HEX(HLA Extractor) is a simple HLA typing tool. It is designed to ???


## Installation and dependencies
It depends on 
- bwa-mem
- Python3
- Perl
- ncbi-blast


## Workflow

### Step 1. Building the database for HLA
`hex_build.py`

By default with the HEX comes the preprocessed database for IMGT <name here.gz> from <date here>.

#### 1a. Loading the IMGT database with HLA nucleotide sequences

1. ftp
2. processing to the fasta file with storing metadata

1. HLA 1st type -> 2-3-4 exons only
2. HLA 2nd type -> 3-4-5 exons only (?)

#### 1b. Construction of the C-value matrix for alleles

C-value ("confuse value") is the probability that one allele is replaced with another by an error.

estimating 1st and 5th exons for missing data


### Step 2. HLA extraction
`hex_extract.py`

#### 2a. Preprocessing input .fastq files

1. filtering out sequences with the bad mean quality
2. simple error correction
3. complex error correction like in BayesHammer? http://www.biomedcentral.com/1471-2164/14/S1/S7

#### 2b. Generate HLA candidates with P-values

1. Cluster and find consensuses
2. Align to the DB
