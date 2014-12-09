function op = rpmst(x,eps)
    m = length(x);
    n = length(eps);
    nbigs = mod(n,m);
    
    p = ceil(n/m);

    
    op = repmat(x',p,1);
    op = reshape(op,m*p,1);
    
    if nbigs == 0 
        op = sqrt(m/n)*op.*eps;
    else
        idx = 1:p*m;
        idx = (mod(idx,p) ~= p-1 ) | (idx./p <= nbigs);
        op = sqrt(m/n)*op(idx).*eps;
    end

       
end
        