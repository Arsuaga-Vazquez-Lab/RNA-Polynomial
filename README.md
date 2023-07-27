# RNA-Polynomial

## Overview

This repository includes code and data for experiments in the following paper. 

[Pengyu Liu, Jacob Lusk, Mariel Vázquez, Analyzing the link between RNA secondary structures and R-loop formation with tree polynomials, *Preprint* (2023).](https://)

## Setup

Running the code requires [Matlab](https://matlab.mathworks.com) with the Bioinformatics Toolbox and the Statistics and Machine Learning Toolbox; see [Matlab toolboxes](https://www.mathworks.com/products.html).

After installing Matlab and required Matlab toolboxes, 
1. clone (or download) this repository to a local directory;
2. go to the folder *data* and unzip the compressed files *pfc8_snrpn.zip* and *pfc53_airn.zip*;
3. keep the unzipped folders *pfc8_snrpn* and *pfc53_airn* in the folder *data*.

Then, the setup is complete.

## Usage

### To compute the tree polynomial representations of an RNA secondary structure

The function *[functions/dbn2pl.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/functions/dbn2pl.m)* converts a dot-bracket notation of an RNA secondary structure without pseudoknots into a parent list of its tree representation of a selected type.
The function *[functions/pl2poly.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/functions/pl2poly.m)* then computes the corresponding polynomial with the parent list of the tree representation.
See the file *[Example.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Example.m)* for an example of using the functions to compute the tree polynomial representations of an RNA secondary structure.

To compute the tree polynomial representations of one's own data of RNA secondary structures,
1. convert an RNA secondary structure to the dot-bracket notation `dbn`;
2. use `L = dbn2pl(dbn,type)` convert the dot-bracket notation `dbn` to the parent list `L` of its tree representation of a selected `type`;
3. use `P = pl2poly(dbn,type)` to compute the corresponding polynomial `P` from the parent list `L` of the tree representation.



### To reconduct the experiments and reproduce the results

For clustering non-coding RNA secondary structures, run the file *[Process_Clustering.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Process_Clustering.m)*.
This will compute the tree polynomial representations of the RNA secondary structures with dot-bracket notations stored in the folder *[data/bpRNA-Rfam](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/tree/main/data/bpRNA-Rfam)* and compute the pairwise polynomial distances between the RNA secondary structures.
The results are stored in the file *[data/bpRNA_Rfam.mat](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/bpRNA_Rfam.mat)*.

For analyzing the link between RNA secondary structures and R-loop formation, run the file *[Process_Rloop.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Process_Rloop.m)*.
This will compute the polynomial representations of the RNA secondary structures produced by DrTransformer ([repository](https://github.com/ViennaRNA/drtransformer), [paper](https://doi.org/10.1093/bioinformatics/btad034)) with dot-bracket notations stored in the folder *data/pfc8_snrpn* and *data/pfc53_airn* and compute the coefficient sums.
The results are stored in the files with names like *data/pFC8_type2.mat*.

If one would like to reproduce the RNA secondary structures with DrTransformer, then one should install DrTransformer (see the instructions in its [repository](https://github.com/ViennaRNA/drtransformer)) and use the *[drtransformer-runner script](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/tree/main/drtransformer-runner)* which segments the sequences of the plasmids stored in the files *[data/pfc8_snrpn_coding_strand.fa](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc53_airn_coding_strand.fa)* and *[pfc53_airn_coding_strand.fa](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc53_airn_coding_strand.fa)* and calls DrTransformer for predicting the secondary structrues of the sequence segments.
This will generate the files of dot-bracket notations in the folder *data/pfc8_snrpn* and *data/pfc53_airn*.

### To view the results in the paper and additional results

For clustering non-coding RNA secondary structures, run the file *[Results_Clustering.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Results_Clustering.m)*.
This will use the file *[data/bpRNA_Rfam.mat](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/bpRNA_Rfam.mat)* to perform tSNE visualization and k-medoids clustering with the pairwise polynomial distances between the non-coding RNA secondary structures and compute the misclassification rates.

For visualizing the correlations between the scaled sums and the probabilities of R-loop formation, run the file *[Results_Correlations.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Results_Correlations.m)*.
This will use the files with names like *data/pFC8_type2.mat* to compute the scaled sums and compare with probabilities of R-loop formation constructed from the files with names like *data/pFC8_SUPERCOILED.bed*.

For analyzing RNA secondary structures that contribute to large coefficient sums, run the file *[Results_Structures.m](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/Results_Structures.m)*.
This will use the files with names like *data/pFC8_type2.mat* to locate the RNA segments with the largest coefficient sums and visualize the RNA secondary structures predicted by DrTransformer with dot-bracket notations stored in the folder *data/pfc8_snrpn* and *data/pfc53_airn*.






