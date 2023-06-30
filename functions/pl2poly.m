function [P,Y] = pl2poly(L,type)
%PL2POLY computes the polynomial representation P for the tree coded by
%the parent list L. The optional output Y records polynomials at each node.

addpath("functions/ndSparse/")

%choose a polynomial based on the tree type
if type == 2 || type == 4

    option = 2; %type 3 and 6 should use polynomial Q

else

    option = 1; % other types should use polynomial P
    
end


switch option

    case 1 %polynomial P

        n = size(L,1); %get the number of nodes

        for i = 1:n %find all tips

            if ~any(L(:,2) == L(i,1))

                L(i,3) = 1; %assign level 1 to tips

            end

        end

        while any(L(:,3) == 0) %assign levels to internal nodes

            for i = 1:n

                if L(i,3) == 0 && ~any(L(L(:,2) == L(i,1),3) == 0)

                    L(i,3) = max(L(L(:,2) == L(i,1),3)) + 1;

                end

            end

        end

        Y = cell(n,1); %preallocate polynomials for all nodes

        for i = 1:n

            Q = zeros(n,n);

            if L(i,3) == 1 %set initial polynomials for tips
                Q(1,2) = 1;
            end

            Y{i} = sparse(Q);

        end

        m = max(L(:,3)); %get maximum level

        for i = 2:m %compute polynomials for each level

            Ni = L(L(:,3)==i,1); %find all nodes at level i

            ni = length(Ni); %count the number of level-i nodes

            for j = 1:ni

                Temp = zeros(n,n); %initialize temorary polynomial
                Temp(1,1) = 1;

                Nj = L(L(:,2)==Ni(j),1); %find child nodes
                nj = length(Nj);

                for k = 1:nj %multiply the polynomials for child nodes

                    Temp = ifftn(fftn(Temp).*fftn(full(Y{Nj(k)})));
                    Temp = round(Temp);

                end

                Temp(2,1) = Temp(2,1) + 1; %plus y

                Y{Ni(j)} = sparse(Temp); %update polynomial

            end

        end

        P = Y{n};

    case 2 %Polynomial Q

        n = length(L); %get the number of nodes

        for i = 1:n %find all tips

            if ~any(L(:,2) == L(i,1))

                L(i,4) = 1; %assign level 1 to tips

            end

        end

        while any(L(:,4) == 0) %assign levels to internal nodes

            for i = 1:n

                if L(i,4) == 0 && ~any(L(L(:,2) == L(i,1),4) == 0)

                    L(i,4) = max(L(L(:,2) == L(i,1),4)) + 1;

                end

            end

        end

        Y = cell(n,1); %preallocate polynomials for all nodes

        for i = 1:n

            Q = zeros(n,n,n);

            if L(i,4) == 1 %set initial polynomials for tips
                Q(1,2,1) = 1;
            end

            Y{i} = ndSparse(Q);

        end

        m = max(L(:,4)); %get maximum level

        for i = 2:m %compute polynomials for each level

            Ni = L(L(:,4)==i,1); %find all nodes at level i

            ni = length(Ni); %count the number of level-i nodes

            for j = 1:ni

                Temp = zeros(n,n,n); %initialize temorary polynomial
                Temp(1,1,1) = 1;

                Nj = L(L(:,2)==Ni(j),1); %find child nodes
                nj = length(Nj);

                for k = 1:nj %multiply the polynomials for child nodes

                    Temp = ifftn(fftn(Temp).*fftn(full(Y{Nj(k)})));
                    Temp = round(Temp);

                end

                if L(L(:,1)== Ni(j),3) == 1 %artificial nodes

                    Temp(1,1,2) = Temp(1,1,2) + 1; %plus z

                elseif L(L(:,1)== Ni(j),3) == 0 % true nodes

                    Temp(2,1,1) = Temp(2,1,1) + 1; %plus y

                end

                Y{Ni(j)} = ndSparse(Temp); %update polynomial

            end

        end

        P = Y{n};

end

end