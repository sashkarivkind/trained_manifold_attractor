dd=2; %direction - 1 - transversal, 2- parallel
Meff=length(pts);
r=sim.r(:,pts);
K=6;
g=net.g;
sim.rp=net.phip(x);
rp=sim.rp(:,pts);
x1full=(x-net.wfb*sim.f_ol);
x1=x1full(:,pts);

Crr=r'*r;
[eta,lambda]=eig(Crr);
lambda=1/Meff*1/net.N*lambda; %doublecheck normalization

[lambda,ii]=sort(diag(lambda),'descend');
lambda=abs(real(lambda))+1e-20;
eta=eta(:,ii);
eta0=diag(eta(:,1));

lambda=lambda(1:K);
eta=eta(:,1:K);

v=r*eta;
v=v*inv(diag(sqrt(lambda)));

B0=1/Meff*1/net.N*eta'*((diag(rp(:,1))*r).'*net.wfb(:,dd));
% B0=1/Meff*1/net.N*eta'*(diag(net.wfb(:,dd))*(r)).'
% B0=1/Meff*1/net.N*eta'*(diag(rp(:,1))*diag(net.wfb(:,dd))*r).'

an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization

B1=1/Meff*1/net.N*eta'*(diag(rp(:,1))*r).'*an; %%doublecheck for transpositions
B1=B1;
zzz=1;
alpha=g^2*inv(zzz*diag(g*sqrt(lambda))-g^2*B1)*B0/zzz;

qn=v.'*net.wout;
Gn=zzz*g^(-1)*sqrt(lambda).*alpha;
G=qn'*Gn

% figure; plot(x1(:,5),an*diag(sqrt(lambda))*eta(5,:)','x')

%% debug:
X1=X-net.wfb(:,2);
lhs=1/net.N*X1'*x1;
rhs=g^2*1/net.N*r'*(rp(:,1).*X);
figure; plot(lhs,rhs,'x')

% %%
% lhs=1/net.N*X1'*(x1*eta);
% rhs=g^2*1/net.N*(rp(:,1).*X)'*(r*eta);
% figure; plot(lhs,rhs,'x')
% 
% %% plugging in a's
% lhs=1/net.N*X1'*(g*an*diag(sqrt(Meff*lambda)));
% rhs=g^2*1/net.N*(rp(:,1).*X)'*(r*eta);
% figure; plot(lhs,rhs,'x')
% 
% %% reopening X into X1+wfb
% lhs=1/net.N*X1'*(g*an*diag(sqrt(Meff*lambda)));
% rhs=g^2*1/net.N*(rp(:,1).*(X1+net.wfb(:,dd)))'*(r*eta);
% figure; plot(lhs,rhs,'x')

%% equation
lhs=(g*diag(sqrt(Meff*lambda)));
rhs0=g^2*1/net.N*(rp(:,1).*(net.wfb(:,dd)))'*(r*eta);
rhs1=g^2*1/net.N*(diag(rp(:,1))*(an))'*(r*eta);
alphaDebu = inv(lhs-rhs1')*rhs0';
%% verifying alphas
lhs=1/net.N*(an*alphaDebu)'*(g*an*diag(sqrt(Meff*lambda)));
rhs=g^2*1/net.N*(rp(:,1).*(an*alphaDebu+net.wfb(:,dd)))'*(r*eta);
% figure; plot(lhs,rhs,'x')
%% recovering gain
% qnDebu=an'*net.wout;
rpc=r*eta;
qnDebu=net.wout'*rpc/diag(lambda)/Meff;
GnDebu=zzz*g^(-1)*sqrt(Meff*lambda).*alphaDebu;
Gdebu=qnDebu*GnDebu