clear;close all;
dd=1; %direction - 1 - transversal, 2- parallel
a=2;b=1;h=0.8;
K=12;
net.g=0.2;
rot=pi/a;
% rot=0;

% hp.alpha_reg=1e-10;
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

%simulating closed loop
x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
r_test_rand=net.phi(x_test_rand);
sim.z_test_rand= net.wout'*r_test_rand;

x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
r_test_noi=net.phi(x_test_noi);
sim.z_test_noi= net.wout'*r_test_noi;

%plotting results
figure;
plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
hold on;
plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);


pt=1;
ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
ee_cl1=eig((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt))));
ee_cl2=eig((net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(x(:,pt))));
% ee_clg0=eig((0*net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));

figure;
plot(ee_cl-1,'*');
hold on;

plot(real(ee_cl1(1:5)-1),imag(ee_cl1(1:5)-1),'o');
plot(real(ee_cl2(1:2)-1),imag(ee_cl2(1:2)-1),'s');

legend( 'full','transversal','tangent');
%

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
lhs=(g*diag(sqrt(Meff*lambda)));
B0=g^2*1/net.N*(rp(:,pt).*(net.wfb(:,dd)))'*(r*eta);
B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*eta);
B0=B0';
B1=B1';
alpha= inv(zzz*lhs-B1)*B0/zzz;
% recovering gain
rpc=r*eta;
qn=net.wout'*rpc/diag(lambda)/Meff;
Gnprefac=g^(-1)*sqrt(Meff*lambda);
Gn=zzz*Gnprefac.*alpha;
MM=inv(lhs)* (B1+B0*((qn(dd,:).*Gnprefac')));
G=qn*Gn
y=eig(MM)-1;
hold on;plot(y,'x','MarkerSize',15)
figure;imagesc(MM)
figure;imagesc(inv(lhs)*B1);title('beta1')



%% linearized net
% lin_net=net;
% lin_net.W=lin_net.W*diag(lin_net.phip(x(:,1)));
% lin_net.phi  = @(x) x;
% lin_net.phip = @(x) ones(size(x));
% X =  fast_conv_to_fp(lin_net,[0;1],struct('ol',1));