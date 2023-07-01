# RNA-Polynomial

## Overview

This repository includes code and data for experiments in the following paper. 

[Pengyu Liu, Jacob Lusk, Mariel Vázquez, Analyzing the link between RNA secondary structures and R-loop formation with tree polynomials, *Preprint* (2023).](https://)

## Setup

Running the code requires [Matlab](https://matlab.mathworks.com) and additional [Matlab toolboxes](https://www.mathworks.com/products.html) including the Bioinformatics Toolbox and the Statistics and Machine Learning Toolbox.

After installing Matlab and required Matlab toolboxes, 
1. clone (or download) this repository to a local directory;
2. go to the folder *data* and unzip the compressed files *pfc8_snrpn.zip* and *pfc53_airn.zip*;
3. keep the unzipped folders *pfc8_snrpn* and *pfc53_airn* in the folder *data*.

Then, the setup is complete.

## Usage

### To compute the tree polynomial representations of an RNA secondary structure

The function *functions/dbn2pl.m* converts a dot-bracket notation of an RNA secondary structure without pseudoknots into a parent list of its tree representation of a selected type.
The function *functions/pl2poly.m* then computes the corresponding polynomial with the parent list of the tree representation.
See the file *Example.m* for an example of using the functions to compute the tree polynomial representations of an RNA secondary structure.

To compute the tree polynomial representations of your own data of an RNA secondary structure,
1. convert the RNA secondary structure to the dot-bracket notation `dbn`;
2. use `L = dbn2pl(dbn,type)` convert the dot-bracket notation `dbn` to the parent list `L` of its tree representation of a selected `type`;
3. use `P = pl2poly(dbn,type)` to compute the corresponding polynomial `P` from the parent list `L` of the tree representation.



### To reconduct the experiments and reproduce the results

For clustering non-coding RNA secondary structures, run the file *Process_Clustering.m*.
This will compute the tree polynomial representations of the RNA secondary structures with dot-bracket notations stored in the folder *data/bpRNA-Rfam* and compute the pairwise polynomial distances between the RNA secondary structures.
The results are stored in the file *data/bpRNA_Rfam.mat*.

For analyzing the link between RNA secondary structures and R-loop formation, run the file *Process_Rloop.m*
This will compute the polynomial representations of the RNA secondary structures produced by DrTransformer ([repository](https://github.com/ViennaRNA/drtransformer), [paper](https://doi.org/10.1093/bioinformatics/btad034)) 


### To view the results in the paper and additional results












