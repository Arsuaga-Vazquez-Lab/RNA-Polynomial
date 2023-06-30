%% Processing data and visualization for R-loops

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

    case 2

        sp = 81; %pFC53 AIRN gene starting point
        ep = 1829; %pFC53 AIRN gene ending point

        pname = 'pFC53';
        plen = 3906;
        pid = 53;

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

BedGY = strcat(pname,'_',tnameGY,'.bed');
FreqGY = bed2freq(BedGY,plen,1,pid);
FreqGY = FreqGY(ep:-1:sp);

%load coefficient sums data
data = strcat(pname,'_type',num2str(type));
load(data)

%scale coefficient sums
ScaledSums = cs2ss(CoeffSums,seqlen,seglen,nseg);

lastten = 1;

[CorSC,StpSC] = mostfit(ScaledSums,FreqSC,seglen,lastten)
[CorGY,StpGY] = mostfit(ScaledSums,FreqGY,seglen,lastten)

%plotting
fs = 16; %set font size
lw = 2; %line width
Colors = lines(7);

figure('Position', [0 1000 1150 650]);

subplot(2,1,1)

hold on
grid on
box on

plot(1:seqlen,FreqSC,'Color',Colors(1,:),'LineWidth',lw);
plot(1:seqlen,ScaledSums(StpSC,:).*max(FreqSC)./max(ScaledSums(StpSC,:)),'Color',[0,0,0],'LineWidth',lw)

xlabel("Position (bp) in the coding strand (5' to 3')")
ylabel('Probability / scaled sum')

xticks(0:200:1800)
ylim([0,0.8])

lgd1 = 'R-loop formation probability (supercoiled)';
lgd2 = strcat(['Scaled sum at transcriptional step',' ',num2str(StpSC)]);
legend(lgd1,lgd2,'Location','northeast')

ttl = strcat([pname,',',' ',ttlSC,',',' ','type',' ',num2str(type),' ','tree representation,',' ',polyname]);
title(ttl)

set(gca,'fontname','Palatino','fontsize',fs)

subplot(2,1,2)

hold on
grid on
box on

plot(1:seqlen,FreqGY,'Color',Colors(2,:),'LineWidth',lw);
plot(1:seqlen,ScaledSums(StpGY,:).*max(FreqGY)./max(ScaledSums(StpGY,:)),'Color',[0,0,0],'LineWidth',lw)

xlabel("Position (bp) in the coding strand (5' to 3')")
ylabel('Probability / scaled sum')

xticks(0:200:1800)
ylim([0,0.8])

lgd1 = 'R-loop formation probability (gyrase)';
lgd2 = strcat(['Scaled sum at transcriptional step',' ',num2str(StpGY)]);
legend(lgd1,lgd2,'Location','northeast')

ttl = strcat([pname,',',' ',ttlGY,',',' ','type',' ',num2str(type),' ','tree representation,',' ',polyname]);
title(ttl)

set(gca,'fontname','Palatino','fontsize',fs)

%% Export data for 3D visualization

writematrix(ScaledSums,strcat(data,'.csv'))
