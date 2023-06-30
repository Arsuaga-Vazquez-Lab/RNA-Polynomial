function [F,Ba] = bed2freq(File,n,nor,id)
%BED2FREQ turns a bed file of r-loops to the vector of frequency of r-loops
%at each nucelotide position.

C = importdata(File);
c = C{1};
if id == 8
    del = c(5);
elseif id == 53
    del = c(6);
end

opts = detectImportOptions(File,"FileType","text");
opts.Delimiter = del;
B = readtable(File,opts);

if nor % if normalize

    B = B(:,2:4);

    B1 = B.Var4;
    NB = zeros(length(B1),1);
    nb = length(B1);

    for i = 1:nb

        sp = strsplit(B1{i},'/');
        NB(i) = str2double(sp{2});

    end

    B.Var4 = NB;

    U1 = unique(NB);
    B.Var5 = zeros(nb,1);

    for i = 1:length(U1)

        ts = min(B.Var2(B.Var4 == U1(i)));

        B.Var5(B.Var2 == ts & B.Var4 == U1(i)) = 1;

    end
    Ba = B;
    Br = B(B.Var5==1,:);



    a = Br.Var2;
    b = Br.Var3;
    m = length(a);

    F = zeros(m,n);

    for i = 1:m

        F(i,a(i):b(i)) = 1;

    end
    F = sum(F);

    F = F./nb;

else

    a = B.Var2;
    b = B.Var3;
    m = length(a);

    F = zeros(m,n);

    for i = 1:m

        F(i,a(i):b(i)) = 1;

    end
    F = sum(F);

    F = F./max(F);

end

end