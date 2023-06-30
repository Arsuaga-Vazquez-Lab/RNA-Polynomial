# RNA-Polynomial

## Overview

This repository includes code and data for experiments in the following paper. 

[Pengyu Liu, Jacob Lusk, Mariel Vázquez, Analyzing the link between RNA secondary structures and R-loop formation with tree polynomials, *Preprint* (2023).](https://)

## Usage

The function **dbn2pl.m** converts a dot-bracket notation of an RNA secondary structure without pseudoknots into a parent list of its tree representation of a selected type.

The function **pl2poly.m** then computes the corresponding polynomial with the parent list of the tree representation.

See the file **Example.m** for an example of using the functions to compute the tree polynomial representations of an RNA secondary structure.
