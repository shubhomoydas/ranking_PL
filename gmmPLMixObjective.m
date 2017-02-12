function [ val ] = gmmPLMixObjective( Pd, G, Pi, D, N, K )
%GMMPLMIXOBJECTIVE Summary of this function goes here
%   Detailed explanation goes here

val = 0;
for k = 1:K
    S  = zeros(N,N);
    sum_pk = sum(Pi(k,:));
    for d = 1:D
        S = S + (Pi(k,d)/sum_pk) .* Pd{d};
    end
    t = S*G(:,k);
    %t'*t
    val = val + t'*t;
end

end

