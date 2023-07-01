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

The function *dbn2pl.m* converts a dot-bracket notation of an RNA secondary structure without pseudoknots into a parent list of its tree representation of a selected type.
The function *pl2poly.m* then computes the corresponding polynomial with the parent list of the tree representation.
See the file *Example.m* for an example of using the functions to compute the tree polynomial representations of an RNA secondary structure.

To compute the tree polynomial representations of your own data of RNA secondary structures,
1. convert the RNA secondary structures to dot-bracket notations;
2. use 'L = dbn2pl(dbn,type)' convert the dot-bracket notation to the parent list of its tree representation of a selected type;
3. use 'P = pl2poly(dbn,type)'to compute the corresponding polynomial of the tree representation.

### To reconduct the experiment and reproduce the results



### To view the results in the paper and additional results










