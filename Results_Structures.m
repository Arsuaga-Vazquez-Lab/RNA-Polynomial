%% Plot secondary structures

clear,clc

addpath("functions/ndSparse/")
addpath("functions/")
addpath("data/")

plasmid = 1; %choose plasmid: 1 - pFC8 and 2 - pFC53

type = 1; %select tree type

switch plasmid

    case 1

        sp = 81; %pFC8 SNRPN gene starting point
        ep = 1512; %pFC8 SNRPN gene ending point

        pname = 'pFC8';
        plen = 3669;
        pid = 8;

        faname = 'pfc8_snrpn_coding_strand.fa';
        foname = 'data/pfc8_snrpn/pfc8_snrpn_coding_strand_';

    case 2

        sp = 81; %pFC53 AIRN gene starting point
        ep = 1829; %pFC53 AIRN gene ending point

        pname = 'pFC53';
        plen = 3906;
        pid = 53;

        faname = 'pfc53_airn_coding_strand.fa';
        foname = 'data/pfc53_airn/pfc53_airn_coding_strand_';

end


tnameSC = 'SUPERCOILED';
tnameGY = 'GYRASE';

ttlSC = 'linearized supercoiled plasmid';
ttlGY = 'linearized gyrase plasmid';


if type == 2 || type == 4

    polyname = 'polynomial Q';

else

    polyname = 'polynomial P';

end

seqlen = ep - sp + 1; %length of gene sequence
seglen = 200; %segment length

nseg = seqlen - seglen + 1; %number of segments


%load bed files and turn into c->t frequency
BedSC = strcat(pname,'_',tnameSC,'.bed');
FreqSC = bed2freq(BedSC,plen,1,pid);
FreqSC = FreqSC(ep:-1:sp);

%load coefficient sums data
data = strcat(pname,'_type',num2str(type));
load(data)

%scale coefficient sums
[ScaledSums,Pos,MA] = cs2ss(CoeffSums,seqlen,seglen,nseg);

lastten = 1;

[CorSC,StpSC] = mostfit(ScaledSums,FreqSC,seglen,lastten);

Seq = fastaread(faname);
Seq = Seq.Sequence;

%% Select segments

%prepare data for segment positions
CS = CoeffSums(StpSC,:);
CS = CS';
CS(:,2:3) = Pos(:,1:2);

%sort coefficient sums
SCS = sortrows(CS,1);

%choose the largest coefficient sums
sel = size(SCS,1)-0:-1:size(SCS,1)-1;

%information about the secondary structurs of the segments
pos = [SCS(sel,2:3),SCS(sel,1)./SCS(end,1)]
StpSC

%record number of secondary structures
nn = length(sel);


%% Plot structures

fs = 16; %set font size
lw = 2; %line width
Colors = lines(7);

for j = 1:nn

    F = strcat(foname,num2str(pos(j,1)-1),'-',num2str(pos(j,2)),'.drf');
    T = readtable(F,'FileType','text','ReadVariableNames',false);

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

    dbn = DBN{StpSC};

    W = Seq(pos(j,1):(pos(j,1)+StpSC-1));

    rnaplot(dbn,'Sequence',strrep(W,'t','u'),'Format','Diagram')
    set(gca,'fontname','Palatino','fontsize',fs)
    set(gcf,'Position',[0 1000 650 500])

    L = dbn2pl(dbn,type);

    nl = size(L,1); %number of nodes

    %generate the adjacency matrix from the parent list
    A = sparse(L(1:end-1,1),L(1:(end-1),2),1);
    A(nl,:) = 0;
    A = A+A';
    G = graph(A);

    %plot the graph object
    figure

    set(gcf,'position',[0,0,700,350]);
    h = plot(G,'LineWidth',2,'Layout','layered','Sources',nl,'MarkerSize',7,'NodeLabel',{});
    h.EdgeColor=[0,0,0];
    h.NodeColor=[0,0,0];

end
