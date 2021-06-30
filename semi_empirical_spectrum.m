function o=semi_empirical_spectrum(net,sim,K,pts,pt,dd,rotate_eigenvectors)

if nargin<6
    dd=[1,2];
end
if nargin<7
    rotate_eigenvectors=0;
end
Meff=length(pts);
x=sim.x;
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

lambda=lambda(1:K);
eta=eta(:,1:K);

an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization
zzz=1;

% equation
lhs=(g*diag(sqrt(Meff*lambda)));
B0=g^2*1/net.N*(diag(rp(:,pt))*(net.wfb(:,:)))'*(r*eta);
B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*eta);
B0=B0';
B1=B1';
alpha= inv(zzz*lhs-B1)*B0/zzz;
% recovering gain
rpc=r*eta;
qn=net.wout'*rpc/diag(lambda)/Meff;
Gnprefac=g^(-1)*sqrt(Meff*lambda);
Gn=zzz*diag(Gnprefac)*alpha;
MM=inv(lhs)* (B1+B0(:,dd)*((qn(dd,:)*diag(Gnprefac))));
G=qn*Gn;
o.G=G;
o.B1=B1;
o.B0=B0;
o.MM=MM;
o.lhs=lhs;
o.Gnprefac=Gnprefac;
o.qn=qn;
o.alpha=alpha;
o.eta=eta;

%% approximating tau eff as the pole that corresponds to the leading eignevalue:

[vv,ee] = eig(G);
[ee,ii]=sort(diag(ee),'descend');
vv = vv(:,ii);
% ttt=eye(size(lhs))-inv(lhs)*B1;
% o.tau_inv_mat=vv\ttt*vv;
% o.tau_inv=o.tau_inv_mat(1,1);

%% use d/ds G(s) @0 = -tau to compute tau
ds=0.01;
zzzds=1+ds;
alphads= inv(zzzds*lhs-B1)*B0/zzzds;
Gnds=zzzds*diag(Gnprefac)*alphads;
Gds=qn*Gnds;
tau_mat_ds = -vv\(Gds-G)*vv/ds;
o.tau_approx_by_deriv=tau_mat_ds(1,1);
0;

