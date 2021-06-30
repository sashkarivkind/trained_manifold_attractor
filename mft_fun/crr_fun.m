function o=crr_fun(sigma,beta_out,psi_vec,A,phi)
o=[];
for ii=1:length(beta_out)
    c = beta_out(ii)/sigma^2;
    b=real(sqrt(1-c^2));
    psi = psi_vec(ii);
    i1 = @(y,t) phi(sigma*(c*y(:,1)+b*y(:,2))+...
    A*(cos(2*pi*t(:,1))*cos(psi)+sin(2*pi*t(:,1))*sin(psi))).*...
    phi(sigma*y(:,1)+A*cos(2*pi*t(:,1)));
    o(end+1) = num_expectation_gauss_uni(i1,2,1,1e6);
end