function o=bulk_fun(g,sigma,A,phip)
o=[];
i1 = @(y,t) phip(sigma*y(:,1)+...
    A*(cos(2*pi*t(:,1)))   ).^2;
o(end+1) = g.*sqrt(num_expectation_gauss_uni(i1,2,1,1e6))-1;
end