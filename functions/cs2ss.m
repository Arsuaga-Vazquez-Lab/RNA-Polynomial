function [SS,Pos,MA] = cs2ss(CS,seqlen,seglen,nseg)
%CS2SS computes the scaled sum from a coeffcient sum.

Pos = zeros(nseg,3); %segment positions
for i = 1:nseg

    Pos(i,1) = i; %segment starting point
    Pos(i,2) = i + seglen - 1; %segment ending point
    Pos(i,3) = (2*i + seglen - 1)/2; %segment center;

end

SS = zeros(seglen,seqlen); %preallocate sclaed sum

for j = 1:seglen

    R = CS(j,:); %get the coefficient sums for transcriptional step j

    A = zeros(seqlen,1); %preallocate overlapping sums for the sequence

    for i = 1:nseg %adding every segment

        ta = zeros(seqlen,1);
        ta(Pos(i,1):Pos(i,2)) = R(i);

        A = A + ta;

    end

    MA = max(A);
    A = A./MA; %scaled to 1

    SS(j,:) = A;

end

end