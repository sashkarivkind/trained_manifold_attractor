%define Francescian connectivity
xi=randn(1000,2);
W=xi*xi';

%extract connectivity parameters
theta=atan2(xi(:,2),xi(:,1));
rho = sqrt(sum(xi.^2,2));
dtheta=repmat(theta,1,1000)-repmat(theta',1000,1);

%distance dependent connectivity plots
figure; plot(dtheta(:),W(:),'.')
[~,uu]=sort(theta);
figure; imagesc(W(uu,uu))
figure; imagesc(inv(diag(rho(uu)))*W(uu,uu)*inv(diag(rho(uu))))