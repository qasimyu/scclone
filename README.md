# SCClone

SCClone is a software for intra-tumor heterogeneity inference from single-cell DNA sequencing data.

## Requirements

* Linux systems.
* CMake3.0+.
* g++.

## Installation

To build binary, do as follows:

```
tar -zxvf SCClone.tar.gz
cd SCClone
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

### Genotype Matrix

The SNVs of single cells are denoted as a genotype matrix. Each row defines the mutation states of a single cell, and each column represents one mutation. Columns are separated by tabs. The genotype matrix is binary.

The entry at position [i,j] should be

* 0 if mutation j is not observed in cell i,
* 1 if mutation j is observed in cell i, or
* 3 if the genotype information is missing

## Output Files

The base name of the output files is provided by users.

### Subclone genotypes

The subclone genotypes are written to a file with suffix "clone_genotypes".

### Cell assignments

The subclone indexs of the cells are written to a file with suffix "cell_assignment".

## Arguments

* `-i, --input <filename>` Replace \<filename\> with the file containing the genotype matrix.

* `-o, --output <string>` Replace \<string\> with the base name of the output file.

## Optional arguments

* `-K, --maxc <INT>` Set \<INT\> to a positive integer. This specifies the maximum number of clones to consider.

* `-a, --alpha <Double>` Set \<Double\> to estimated false positive rate of the single-cell sequencing experiment.

* `-b, --beta <Double>` Set \<Double\> to estimated false negative rate of the single-cell sequencing experiment.

## Contact

If you have any questions, please contact zhyu@nxu.edu.cn.