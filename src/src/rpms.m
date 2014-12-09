function op = rpms(x,m,eps)
    n = length(x);
    sizes = zeros(m,1);
    nbigs = mod(n,m);

    sizes(1:nbigs) = ceil(n/m);
    sizes(nbigs+1:m) = floor(n/m);
    ends = cumsum(sizes);
    
    x = x.*eps;
    x = cumsum(x);
    x = [0;x(ends)];
    x = diff(x);

    op = sqrt(m/n)*x;
    
end
        
