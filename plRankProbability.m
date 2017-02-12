function [ prob ] = plRankProbability( R, G, N )
%PLRANKPROBABILITY Returns the probabilty of ranking under PL dist.

prob = 1;
for i=1:(N-1)
    prob = prob * (G(R(i))/sum(G(R(i:N))));
end

end

