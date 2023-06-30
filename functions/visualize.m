function visualize(L)
%VISUALIZE plots a bianry tree T with tip labels.

n = size(L,1); %number of nodes

%generate the adjacency matrix from the parent list
A = sparse(L(1:end-1,1),L(1:(end-1),2),1); 
A(n,:) = 0;
A = A+A';
G = graph(A);

%plot the graph object
figure

set(gcf,'position',[0,0,2000,2000]);
h = plot(G,'LineWidth',2,'Layout','layered','Sources',n,'MarkerSize',7);
h.EdgeColor=[0,0,0];
h.NodeColor=[0,0,0];


end