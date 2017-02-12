function [ Pairwise ] = convertToPairwise( Ranking )
%CONVERTTOPAIRWISE Converts input row vector of rankings to pairwise
%   square matrix where diagonal elements are -(sum of column)

n = size(Ranking,2);
[~,RI] = sort(Ranking);
Pairwise = zeros(n);
for i=1:n
    r = RI(i);
    for j=1:n
        Pairwise(i,j) = (i ~= j && r < RI(j));
    end
end
for i=1:n
    Pairwise(i,i) = -sum(Pairwise(:,i));
end

end

