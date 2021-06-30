dd=2; %direction - 1 - transversal, 2- parallel

Meff=length(pts);
r=sim.r(:,pts);
K=12;
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

an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization
zzz=1;

%% equation
lhs=(g*diag(sqrt(Meff*lambda)));
B0=g^2*1/net.N*(rp(:,pt).*(net.wfb(:,dd)))'*(r*eta);
B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*eta);
B0=B0';
B1=B1';
alpha= inv(zzz*lhs-B1)*B0/zzz;
%% recovering gain
rpc=r*eta;
qn=net.wout'*rpc/diag(lambda)/Meff;
Gnprefac=g^(-1)*sqrt(Meff*lambda);
Gn=zzz*Gnprefac.*alpha;
MM=inv(lhs)* (B1+B0*((qn(dd,:).*Gnprefac')));
G=qn*Gn
eig(MM)-1
