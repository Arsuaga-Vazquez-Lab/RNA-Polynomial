%% Visualization

clear,clc

addpath("functions/ndSparse/")
addpath("functions/")
addpath("data/")

load("data/bpRNA_Rfam.mat")

fs = 16; %set font size
ms = 20; %marker size
Colors = lines(7);

type = 1;

D = Dists{type};
C = cmdscale(D,3);
T = tsne(D);

%set up titles for different types
switch type
    case 1
        t1 = strcat('MDS plot, type 1 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 1 tree, polynomial P distance');
    case 2
        t1 = strcat('MDS plot, type 2 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 2 tree, polynomial P distance');
    case 3
        t1 = strcat('MDS plot, type 3 tree, polynomial Q distance');
        t2 = strcat('tSNE plot, type 3 tree, polynomial Q distance');
    case 4
        t1 = strcat('MDS plot, type 4 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 4 tree, polynomial P distance');
    case 5
        t1 = strcat('MDS plot, type 5 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 5 tree, polynomial P distance');
    case 6
        t1 = strcat('MDS plot, type 6 tree, polynomial Q distance');
        t2 = strcat('tSNE plot, type 6 tree, polynomial Q distance');
    case 7
        t1 = strcat('MDS plot, type 7 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 7 tree, polynomial P distance');
    case 8
        t1 = strcat('MDS plot, type 8 tree, polynomial P distance');
        t2 = strcat('tSNE plot, type 8 tree, polynomial P distance');
end

%plotting
figure('Position', [0 1000 1200 450])

subplot(1,2,1)

box on
grid on
hold on

for i = 1:7

    plot(C(Data(Data(:,2)==i,3),1),C(Data(Data(:,2)==i,3),2),'.','MarkerSize',ms)

end

xlabel('Dim. 1')
ylabel('Dim. 2')

title(t1)

set(gca,'fontname','Palatino','fontsize',fs)


subplot(1,2,2)

box on
grid on
hold on

for i = 1:7

    plot(T(Data(Data(:,2)==i,3),1),T(Data(Data(:,2)==i,3),2),'.','MarkerSize',ms)

end

xlabel('Dim. 1')
ylabel('Dim. 2')

title(t2)

legend('5.8S ribosomal RNA','U1 spliceosomal RNA','U2 spliceosomal RNA', ...
    'Vault RNA','U12 minor spliceosomal RNA','U3 small nucleolar RNA', ...
    '6S/SsrS RNA','Location','northeast')

set(gca,'fontname','Palatino','fontsize',fs)

%% Computing misclassification rates

clear,clc

addpath("functions/ndSparse/")
addpath("functions/")
addpath("data/")

load("data/bpRNA_Rfam.mat")

rt = 1000; %repeating times

MR = zeros(rt,8);

for type = 1:8

    D = Dists{type};

    for i = 1:rt

        i

        cls = kmedoids(D,7);

        cmp = [Data(:,2),cls];
        cmp(:,3) = 0;

        %the majority rule
        for j = 1:7

            MD = mode(cmp(cmp(:,1) == j,2));

            cmp(cmp(:,1)==j,3) = cmp(cmp(:,1) == j,2) == MD;

        end

        %record misclassification rates
        MR(i,type) = 1 - sum(cmp(:,3)) / n;

    end

end

mean(MR)

