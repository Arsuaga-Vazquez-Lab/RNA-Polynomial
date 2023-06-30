function L = dbn2pl(B,type)
%DBN2PL converts a dot-bracket notation string of an RNA secondary
%structure to a parent list of six types of trees.

switch type

    case 1

        B = shrink(B); %make stems of length one

        B = erase(B,'.'); %remove all leaves
        B = replace(B,'()','.'); %identify new leaves

        L = convert_p(B);


    case 2

        B = shrink(B); %make stems of length one

        L = convert_q(B);

    case 3

        B = erase(B,'.'); %remove all leaves
        B = replace(B,'()','.'); %identify new leaves

        L = convert_p(B);

    case 4

        L = convert_q(B);

    case 5

        B = shrink(B); %make stems of length one

        L = convert_p(B);

    case 6

        L = convert_p(B);

    case 7

        B = shrink(B); %make stems of length one

        L = convert_q(B);

        L = prune(L,B);

    case 8

        L = convert_q(B);

        L = prune(L,B);

    otherwise

        error('Tree type should be an integer between 1 and 8.')
end

end

%% Supporting functions

function N = prune(L,B)

l = count(B,'.'); %number of tips

L(L(:,1) <= l,:) = []; %remove current leaves
L = [L(L(:,3)==1,:);L(L(:,3) == 0,:)]; %move new leaves to the top of the list

for i = 1:size(L,1)
    L(L(:,1:2) == L(i,1)) = i; %relabeling the nodes
end

N = L;

end



function C = shrink(B)

C = B;

flag = 0;

while flag == 0

    %find basepairs in the original dot-bracket notation
    [x,y] = find(rnaconvert(C)==1);
    BP = sortrows([x,y]);

    %in case there is no basepairs
    if size(BP,1) <= 1

        break

    end

    %check every basepair
    for i = 1:(size(BP,1)-1)

        tflag = 0;

        %remove basepairs if the stem has more than one basepair
        if BP(i+1,1) == BP(i,1) + 1 && BP(i+1,2) == BP(i,2) - 1

            C(BP(i,:)) = [];
            break

        else

            tflag = 1;

        end

    end

    if tflag == 1

        flag = 1;

    end

end

end



function L = convert_q(B)

l = count(B,'.'); %number of tips
k = count(B,'('); %number of internal nodes

C = strrep(B,')','('); %get tip groups
grp = strsplit(C,'(');

sg = cellfun(@length, grp); %size of each group
sg = sg(sg ~= 0); %remove empty groups

g = length(sg); %number of groups

n = l + k + g + 1; %number of nodes: tips...stem nodes...group nodes...root

L = zeros(n,3); %preallocate parent list indx|parent|isartificial
L(:,1) = (1:n);

%assing parent nodes
ind = 1;

for i = 1:g

    L(ind:(ind+sg(i)-1),2) = l+k+i; %assign parent
    ind = ind + sg(i);

end

%build mapping between tree node index and dot-bracket positions for tips
M(:,1) = (1:l)';
M(:,2) = find(B=='.');


%base pair information from dot-bracket notation
[x,y] = find(rnaconvert(B)==1); %base piar positions
BP = sortrows([x,y]);
BP(:,3) = ((l+1):(l+k))'; %tree node correspondence

%build key sequence
m = length(B);
seq = zeros(m,2);

%convert dot-bracket notation to key sequence
for i = 1:m

    if B(i) == '.'
        seq(i,1) = M(M(:,2) == i,1);
    elseif B(i) == '('
        seq(i,1) = BP(BP(:,1) == i,3);
    elseif B(i) == ')'
        seq(i,1) = BP(BP(:,2) == i,3);
    end

end

%assign parent nodes to stem and group nodes

cp = n;

for i = 1:m

    if B(i) == '('

        seq(i,2) = cp(end);
        cp = [cp;seq(i,1)];

    elseif B(i) == ')'

        cp(end) = [];
        seq(i,2) = cp(end);

    elseif B(i) == '.'

        seq(i,1) = L(L(:,1)==seq(i,1),2);
        seq(i,2) = cp(end);

    end

end

useq = unique(seq,'rows');

L((l+1):(l+k+g),2) = useq(:,2);
L((l + k +1):end,3) = 1;
L(L(:,2) == 0,3) = 0;

end



function L = convert_p(B)

l = count(B,'.'); %number of tips
k = count(B,'('); %number of internal nodes

if isempty(B)

    L = [1, 0];

else

    n = l + k + 1; %number of nodes: tips...stem nodes...root

    L = zeros(n,2); %preallocate parent list indx|parent
    L(:,1) = (1:n);

    %build mapping between tree node index and dot-bracket positions for tips
    M(:,1) = (1:l)';
    M(:,2) = find(B=='.');


    %base pair information from dot-bracket notation
    [x,y] = find(rnaconvert(B)==1); %base piar positions
    BP = sortrows([x,y]);
    BP(:,3) = ((l+1):(l+k))'; %tree node correspondence

    %build key sequence
    m = length(B);
    seq = zeros(m,2);

    %convert dot-bracket notation to key sequence
    for i = 1:m

        if B(i) == '.'
            seq(i,1) = M(M(:,2) == i,1);
        elseif B(i) == '('
            seq(i,1) = BP(BP(:,1) == i,3);
        elseif B(i) == ')'
            seq(i,1) = BP(BP(:,2) == i,3);
        end

    end

    %assign parent nodes to stem and group nodes

    cp = n;

    for i = 1:m

        if B(i) == '('

            seq(i,2) = cp(end);
            cp = [cp;seq(i,1)];

        elseif B(i) == ')'

            cp(end) = [];
            seq(i,2) = cp(end);

        elseif B(i) == '.'

            seq(i,2) = cp(end);

        end

    end

    useq = unique(seq,'rows');
    L(1:end-1,2) = useq(:,2);

end

end
