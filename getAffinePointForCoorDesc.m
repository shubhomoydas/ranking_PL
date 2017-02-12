function [ Z ] = getAffinePointForCoorDesc( Z_x, Z_y, alpha, D, K )
%GETAFFINEPOINTFORCOORDESC returns alpha*x + (1-alpha)*y

Z = zeros(K,D);
for d = 1:D
    Z(:,d) = alpha*Z_x(:,d) + (1-alpha)*Z_y(:,d);
end

end

