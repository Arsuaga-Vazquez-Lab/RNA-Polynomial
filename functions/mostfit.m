function [M,step] = mostfit(SS,freq,seglen,option)
%MOSTFIT finds the most correlated transcriptional step to the scaled sum.

if option == 1 %the last 10 steps

    st = 191;

elseif option == 0 %all steps

    st = 1;

end

PC = zeros(seglen,1); %Pearson's correlation coefficients

for i = st:seglen

     pc = corrcoef(SS(i,:),freq);

     PC(i) = pc(1,2);

end

M = max(abs(PC));

step = find(PC == M);

if isempty(step)

    step = find(PC == -M);
    M = -M;

end

end