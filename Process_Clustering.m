%% Clustering ncRNA secondary structures from Rfam database

clear,clc

addpath("functions/ndSparse/")
addpath("functions/")
addpath("data/")

% DATA INFORMATION

% 5.8S: 713-773; max: 153. <- 1 1:61
% U1 spliceosomal: 774-873; max: 133. <- 2 62:161
% U2 spliceosomal: 874-1081; max: 170. <- 3 162:369

% Vault: 2036-2108; max: 126. <- 4 370:442
% U12 minor spliceosomal: 2109-2170; max: 179; <- 5 443:504

% U3 small nucleolar: 2941-3027; max: 173. <- 6 505:591 (excluding 2946,2951,2965,2970,2998 with pseudoknots)
% 6S/SsrS: 3028-3176; max: 154. <- 7 592:740

%% Data preparation

Data = [713:1081,2036:2170,2941:3176]';
Data = setdiff(Data,[2946;2951;2965;2970;2998]);
n = length(Data);

%assign data classes
for i = 1:n

    if Data(i,1) <= 773

        Data(i,2) = 1;

    elseif Data(i,1) <= 873

        Data(i,2) = 2;

    elseif Data(i,1) <= 1081

        Data(i,2) = 3;

    elseif Data(i,1) <= 2108

        Data(i,2) = 4;

    elseif Data(i,1) <= 2170

        Data(i,2) = 5;

    elseif Data(i,1) <= 3027

        Data(i,2) = 6;

    elseif Data(i,1) <= 3176

        Data(i,2) = 7;

    end

end

Data(:,3) = (1:n)'; %assign new indices

%% Compute polynomials and distances

Polys = cell(n,8);

for type = 1:8 %type of tree representation

    parfor i = 1:n

        i

        F = strcat('data/bpRNA-Rfam/bpRNA_RFAM_',num2str(Data(i,1)),'.dbn');
        B = readdbn(F);

        P = dbn2pl(B,type);
        Y = pl2poly(P,type);

        Polys{i,type} = Y;

    end

end

Dists = cell(8,1);

for type = 1:8

    Dists{type} = candist(Polys(:,type),type);

end

save('bpRNA_Rfam')

