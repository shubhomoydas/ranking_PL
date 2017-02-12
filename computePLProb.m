function [ prob ] = computePLProb( P, V )
%COMPUTEPLPROB Summary of this function goes here
%   Detailed explanation goes here
    prob = 1;
    n = length(P);
    tP = 1:n;
    tV = V;
    for i=1:(n-1)
        %tP(i:end)
        %tV(i:end)
        t = P(i);
        v = V(t);
        prob = prob * v / sum(tV(i:end));
        tt = find(tP == t);
        tP(tt) = tP(i);
        tP(i) = t;
        tV(tt) = tV(i);
    end;
end

