%% R-loop formation

clear,clc

addpath("functions/ndSparse/")
addpath("functions/")
addpath("data/")

%initial settings
plasmid = 1; %choose plasmid 1 - pFC8, 2 - pFC53
type = 1; %choose tree type from 1 - 8

seglen = 200;

switch plasmid

    case 1

        Data = (1:1:1233)'; %pFC8
        n = size(Data,1); %number of segments

        filename = strcat('pFC8_type',num2str(type));


    case 2

        Data = (1:1:1550)'; %pFC53
        n = size(Data,1); %number of segments

        filename = strcat('pFC53_type',num2str(type));

end

CoeffSums = zeros(seglen,n); %preallocate coefficient sums

parfor j = 1:n

    j

    %read drtransformer files for each segment
    switch plasmid

        case 1

            F = strcat('data/pfc8_snrpn/pfc8_snrpn_coding_strand_',num2str(Data(j)-1),'-',num2str(Data(j)+seglen-1),'.drf');
            T = readtable(F,'FileType','text','ReadVariableNames',false);

        case 2

            F = strcat('data/pfc53_airn/pfc53_airn_coding_strand_',num2str(Data(j)-1),'-',num2str(Data(j)+seglen-1),'.drf');
            T = readtable(F,'FileType','text','ReadVariableNames',false);

    end

    %making a list of steps of structures
    cur = '';
    DBN = cell(seglen,1);
    Len = zeros(seglen,1);

    for i = 1:seglen

        if ~ismember(i,T.Var4)

            cur = strcat(cur,'.');
            DBN{i} = cur;
            Len(i) = length(cur);

        else

            m = min(T.Var3(T.Var4 == i));
            cur = T.Var2(T.Var3 == m & T.Var4 == i);
            cur = cur{1};
            DBN{i} = cur;
            Len(i) = length(cur);

        end

    end

    %compute the coefficent sums
    for i = 1:seglen

        B = DBN{i};

        L = dbn2pl(B,type);
        P = pl2poly(L,type);

        CoeffSums(i,j) = sum(sum(sum(full(P))));

    end

end

save(filename)
