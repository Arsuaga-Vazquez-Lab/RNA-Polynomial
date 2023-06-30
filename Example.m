%% Example
clear,clc

addpath("data/")
addpath("functions/")
addpath("functions/ndSparse/")

%dot-bracket notation of the example
B = '.(((.((.((...))...((((....)).)).)).)))..';

type = 1; %select tree type from 1 to 8.

%convert dot-bracket notation to parent list
L = dbn2pl(B,type);
%plot the tree
visualize(L)

%compute the polynomial
P = pl2poly(L,type);


