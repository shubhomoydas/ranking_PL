function [ obj ] = surrObj( Z, Phi )
%SURROBJ Summary of this function goes here
%   Detailed explanation goes here
    obj = 0;
    K = size(Z,1);
    for k=1:K
        %  This is concave
        %{
        obj(j) = obj(j) + ...
            sqrt(Phi(:,k)'*Phi(:,k)) ...
            * (1/sqrt((Z(k,:)*Z(k,:)'))) ...
            * Z(k,:)*Phi(:,k);
        %}
        % This appears convex
        %%{
        obj = obj + ...
            sqrt(Phi(:,k)'*Phi(:,k)) * sqrt((Z(k,:)*Z(k,:)')) ...
            * (1/(sum(Z(k,:)))^2) ...
            * Z(k,:)*Phi(:,k);
        %%}
        % This is neither convex nor concave
        %{
        for m=1:D
            for n=1:D
                obj(j) = obj(j) + Z(k,m)*Z(k,n)*Phi(m,k)*Phi(n,k);
            end
        end
        obj(j) = obj(j) * (1/(sum(Z(k,:)))^2);
        %}
    end

end

