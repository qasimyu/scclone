# SCClone

SCClone is a software for intra-tumor heterogeneity inference from single-cell DNA sequencing data.

## Requirements

* Linux systems.
* CMake3.0+.
* g++.

## Installation

To build binary, do as follows:

```
tar -zxvf scclone.tar.gz
cd scclone
cmake .
make
```

After the installation, the main program of SCClone is generated in “bin” directory. Type following command if you want to add SCClone to system path:
```
make install
```

## Usage

SCClone uses the genotype matrix derived from single-cell DNA sequencing data to infer clonal architecture.

Example:

```
scclone -i testdata/example.txt -o testdata/example
```

## Input Files

### 1. Genotype Matrix

The SNVs of single cells are denoted as a genotype matrix. Each row defines the mutation states of a single cell, and each column represents one mutation. Columns are separated by tabs. The genotype matrix can be binary/ternary.

#### a) Binary: Only presence/absence of a mutation is distinguished.
The entry at position [i,j] should be

* 0 if mutation j is not observed in cell i,
* 1 if mutation j is observed in cell i, or
* 3 if the genotype information is missing

#### b) Ternary: Heterozygous and homozygous mutations are distinguished.
The entry at position [i,j] should be

* 0 if mutation j is not observed in cell i,
* 1 if heterozygous mutation j is observed in cell i
* 2 if homozygous mutation j is observed in cell i
* 3 if the genotype information is missing

### 2. Cell names (optional)

A file listing the names of the single cells. Each row specifies the name of a single cell.
If no such file is specified, the cells are numbered from 1 to N (N is the number of cells).

### 3. Mutation names (optional)

A file listing the names of the mutations. Each row specifies the name of a mutation.
If no such file is specified, the mutations are numbered from 1 to M (M is the number of mutations).

## Output Files

The base name of the output files is provided by users.

### Subclone genotypes

The subclone genotypes are written to a file with suffix "clone_genotypes".

### Subclonal tree

The subclonal tree is written to a file in GraphViz format.

## Arguments

* `-i, --input <filename>` Replace \<filename\> with the file containing the genotype matrix.

* `-o, --output <string>` Replace \<string\> with the base name of the output file.

## Optional arguments

* `-c, --clabel <filename>` Replace \<filename\> with the path to the file containing the names of the cells.

* `-m, --mlabel <filename>` Replace \<filename\> with the path to the file containing the names of the mutations.

* `-K, --maxc <INT>` Set \<INT\> to a positive integer. This specifies the maximum number of clones to consider.

* `-t, --threads <INT>`  Set \<INT\> to the number of threads to use. The default value is 1.

* `-a, --alpha <Double>` Set \<Double\> to estimated false positive rate of the single-cell sequencing experiment.

* `-b, --beta <Double>` Set \<Double\> to estimated false negative rate of the single-cell sequencing experiment.

* `-A, --max_alpha <Double>`  Set \<Double\> to the desired maximum false positive rate (only effective when learning FPR from data). The default value is 0.05.

* `-B, --max_beta <Double>`  Set \<Double\> to the desired maximum false negative rate (only effective when learning FNR from data). The default value is 0.5.

## Contact

If you have any questions, please contact zhyu@nxu.edu.cn.