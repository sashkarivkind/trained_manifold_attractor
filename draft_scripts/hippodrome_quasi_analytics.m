dd=2; %direction - 1 - transversal, 2- parallel
Meff=length(pts);
r=sim.r(:,pts);

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

B0=1/Meff*1/net.N*eta'*((diag(rp(:,1))*r).'*net.wfb(:,dd));

an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization

B1=1/Meff*1/net.N*eta'*(diag(rp(:,1))*r).'*an; %%doublecheck for transpositions
B1=B1';
zzz=1;
alpha=g^2*inv(zzz*diag(g*sqrt(lambda))-g^2*B1)*B0/zzz;

v=r*eta;
v=v*inv(diag(sqrt(lambda)));
qn=v.'*net.wout;
Gn=zzz*g*sqrt(lambda).*alpha;
G=qn'*Gn