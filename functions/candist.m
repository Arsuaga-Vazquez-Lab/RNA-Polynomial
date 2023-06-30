function D = candist(Polys,type)
%CANDIST computes the Canberra distances between the polyomials and returns
%the distance matrix D.

addpath("functions/ndSparse/")

%choose a polynomial based on the tree type
if type == 2 || type == 4

    option = 2; %type 3 and 6 should use polynomial Q

else

    option = 1; % other types should use polynomial P
    
end

switch option

    case 1

        n = size(Polys,1); %number of polynomials

        S = cell2mat(cellfun(@size,Polys,'UniformOutput',false)); % polynomials sizes

        s = max(max(S)); %max size

        D = zeros(n,n);

        for i = 1:n

            i

            for j = (i+1):n

                P = zeros(s,s);
                Q = zeros(s,s);
                P(1:S(i,1),1:S(i,2)) = Polys{i};
                Q(1:S(j,1),1:S(j,2)) = Polys{j};

                Temp = abs(P-Q)./(P+Q);
                Temp(isnan(Temp)) = 0;
                D(i,j) = sum(Temp,'all');

            end

        end

        D = D + D';

    case 2

        n = size(Polys,1); %number of polynomials

        S = cell2mat(cellfun(@size,Polys,'UniformOutput',false)); % polynomials sizes

        s = max(max(S)); %max size

        D = zeros(n,n);

        for i = 1:n

            p1 = Polys{i};

            parfor j = (i+1):n

                [i j]

                ts = S;

                q1 = Polys{j};

                P = zeros(s,s,s);
                Q = zeros(s,s,s);
                P(1:ts(i,1),1:ts(i,2),1:ts(i,3)) = full(p1);
                Q(1:ts(j,1),1:ts(j,2),1:ts(j,3)) = full(q1);

                Temp = abs(P-Q)./(P+Q);
                Temp(isnan(Temp)) = 0;
                D(i,j) = sum(Temp,'all');

            end

        end

        D = D + D';

end

end