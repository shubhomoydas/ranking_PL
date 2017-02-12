function [ S ] = sampleFromPLDistribution( Gamma )
%SAMPLEFROMPLDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here
    
    N = size(Gamma,1);
    Items = 1:N;
    G = Gamma';
    S = zeros(1,N);
    for i = 1:N
        G = G / sum(G);
        n = find(mnrnd(1,G) > 0);
        S(i) = Items(n);
        G(n) = []; Items(n) = [];
    end

end

