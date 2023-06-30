function [B,S,pk] = readdbn(F)
%READDBN reads .dbn file located as F and returns the dot-bracket notation
%of the RNA secondary structure.

T = readtable(F,'FileType','text','ReadVariableNames',false); %read dbn file
T = T.Var1;
S = T{3}; %sequence
B = T{4}; %dot-bracket notation

pk = ismember('[',B); %check if pseudoknots exist

end