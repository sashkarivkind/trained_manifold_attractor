function psi_hat=psi_decoded(z,f,psi_vec,method)
%Z - 2xN samples
%F - 2xM samples
%PSI_VEC - 1XM values
%METHOD - METHOD for estimating psi_hat
%OUT:
%PSI_HAT 1xN: estimate of psi of Z

if nargin<4
    method = 'nn';
end
v1=ones(1,size(f,2));
psi_hat =zeros(1,size(z,2)); 
for zz=1:size(z,2)
    this_z=z(:,zz);
    if strcmp(method,'nn')
        [~,ii]=min(sum((this_z*v1-f).^2,1));
        psi_hat(zz)=psi_vec(ii);
    else
        error
    end
end
    