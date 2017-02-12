function [ ll ] = computePLLogLik( S, Gs, Pi, Z, alpha_0, beta_0 )
%COMPUTEPLLOGLIK Summary of this function goes here
%   Detailed explanation goes here

    ll = 0;
    N = size(Gs,1);
    K = size(Gs,2);
    D = size(Z,1);
    Gs_sum = zeros(D,N);
    for k=1:K
        for d = 1:D
            Gs_sum(d,:) = cumsum(Gs(S(d,N:-1:1),k));
            Gs_sum(d,:) = Gs_sum(d,N:-1:1);
        end
        for d=1:D
            ll = ll + Z(d,k) * log(Pi(k));
            for i=1:N
                ll = ll + Z(d,k) * log(Gs(S(d,i),k)/Gs_sum(d,i));
            end
        end
        ll = ll + (alpha_0 - 1) * sum(log(Gs(:,k))) - beta_0 * sum(Gs(:,k));
    end

end

