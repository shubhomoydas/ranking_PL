function [ c,ceq,DC,DCeq ] = confungrad( X )
%CONFUNGRAD Summary of this function goes here
%   Detailed explanation goes here

c = -X;
ceq = [];
if nargout > 2
    DC = -eye(size(X,2));
    DCeq = [];
end

end

