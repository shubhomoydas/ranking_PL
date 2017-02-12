function [ G, Z, P ] = getAffinePoint( R, G_x, P_x, G_y, P_y, alpha, D, N, K )
%GETAFFINEPOINT returns alpha*x + (1-alpha)*y

G = zeros(N,K);
for k = 1:K
    G(:,k) = alpha*G_x(:,k) + (1-alpha)*G_y(:,k);
end
P = alpha*P_x + (1-alpha)*P_y;

Z = zeros(K,D);
for d = 1:D
    %Z(:,d) = alpha*Z_x(:,d) + (1-alpha)*Z_y(:,d);
    Z(:,d) = emProbPLMembership( R(d,:), G, P, N, K );
end

end

