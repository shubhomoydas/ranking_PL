function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    y = sum(r,2);
    y(find(y == 0)) = 1;
    r = r ./ repmat(y,1,p);
    %r = r ./ repmat(sum(r,2),1,p);
