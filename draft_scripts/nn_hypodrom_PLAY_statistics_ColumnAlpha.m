clear;close all;
beta_rescaled_rec={};
for jj=1:20

dd=2; %direction - 1 - transversal, 2- parallel
a=2;b=1;h=0.;
K=14;
net.g=1.;
% rot=pi/a;
rot=0;


net.N=1000;

net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
[hp, net, sim] = prep_network_param(struct,net,struct);
% sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];
% rot=pi/a;

rotMat=[cos(rot), -sin(rot); sin(rot), cos(rot) ];
sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];
sim.f_ol=rotMat*sim.f_ol;figure(11);plot(sim.f_ol(1,:),sim.f_ol(2,:),'o'); hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
hp.M=24;

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%learning output weights
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

pt=1;

%%
pts=1:1:length(sim.f_ol);
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

lambda=lambda(1:K);
eta=eta(:,1:K);

an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization
zzz=1;

% equation

tmp=(g*diag(sqrt(Meff*lambda)));
lhs=tmp^2;
B0=g^2*1/net.N*(rp(:,pt).*(net.wfb(:,dd)))'*(r*eta);
B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*eta);
B0=B0';
B1=B1'*tmp;
alpha= inv(zzz*lhs-B1)*B0/zzz;
% recovering gain
rpc=r*eta;
qn=net.wout'*rpc/diag(lambda)/Meff;
% qn_alt=net.wout'*rpc(rpc*rpc');

Gnprefac=(g^(-1)*sqrt(Meff*lambda')*tmp)';
Gn=zzz*Gnprefac.*alpha;
MM=inv(lhs)* (B1+B0*((qn(dd,:).*Gnprefac')));
G=qn*Gn
y=eig(MM)-1;
%%
jj
beta_rescaled_rec{end+1}=inv(lhs)*B1
end

%%
zz=[];
for qq=1:20
    zz(qq,:,:)=abs(beta_rescaled_rec{qq});
end
mean_mat=reshape(mean(zz,1),K,K);
std_mat=reshape(std(zz,[],1),K,K)
figure;
subplot(2,2,1);
imagesc(mean_mat)
title('mean')
colorbar
subplot(2,2,2);
imagesc(std_mat)
title('std')
colorbar

subplot(2,2,3);
imagesc(mean_mat-diag(diag(mean_mat)))
title('mean no diag')
colorbar

subplot(2,2,4);
imagesc(std_mat-diag(diag(std_mat)))
title('std no diag')
colorbar

