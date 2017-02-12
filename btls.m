function [ X1, k, Fx ] = btls( X0, fn, dfn, beta, gamma )
%BTLS Summary of this function goes here
%   Detailed explanation goes here

tol = 10e-4;
X1 = X0;
k = 0;
Fx = [];
while k < 20
    k = k + 1;
    fx = fn(X1);
    dx = dfn(X1);
    dx2 = sum(sum(dx.^2));
    Fx = [Fx;fx];
    if dx2 <= tol
        break
    end
    alpha = 1;
    while fn(X1-alpha*dx) > fx - gamma*alpha*(dx2)
        alpha = alpha*beta;
    end
    X1 = X1-alpha*dx;
    [k fx alpha]
end
fprintf('Found objective in %d steps\n',k);
end

