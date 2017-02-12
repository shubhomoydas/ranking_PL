function [ Z ] = emProbPLMembership( R, G, P, N, K )
%EMPROBMEMBERSHIP Summary of this function goes here
%   Detailed explanation goes here

Z = zeros(K,1);
for k = 1:K
    Z(k) = P(k) * plRankProbability(R, G(:,k), N);
end
Z = Z / sum(Z); % normalize

end

