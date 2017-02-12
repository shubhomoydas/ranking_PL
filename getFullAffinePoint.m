function [ G, Z ] = getFullAffinePoint( G_x, Z_x, G_y, Z_y, alpha, D, N, K )
%GETFULLAFFINEPOINT returns alpha*x + (1-alpha)*y

G = zeros(N,K);
for k = 1:K
    G(:,k) = alpha*G_x(:,k) + (1-alpha)*G_y(:,k);
end

Z = zeros(K,D);
for d = 1:D
    Z(:,d) = alpha*Z_x(:,d) + (1-alpha)*Z_y(:,d);
end

end

