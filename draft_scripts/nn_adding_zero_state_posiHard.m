clear;
result=struct;
tmax=500;
dt=0.25;

net.g=1.5;
net.N=1000;
hp.bb=-0.1;
hp.A=2;
hp.M=40;
phi_tmp=@(x) min(max(x,zeros(size(x))),ones(size(x))); 
phip_tmp =@(x) x>0&x<1;
net.phi=@(x) phi_tmp(x+hp.bb); 
net.phip =@(x) phip_tmp(x+hp.bb);

[hp, net, sim] = prep_network_param(hp, net,struct);
net.J0 = 0;
net.W=net.W+net.J0/net.N;

%finding realistic initial conditions
% x0a = fast_conv_to_fp(net,[1;0]*[0:0.1:3],struct('ol',1));
% x0 = fast_conv_to_fp(net,[1;0]*hp.A,struct('ol',1));

%performing open loop simulation
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1 ...
    ));
%     ,'xinit',repmat(x0,1,size(sim.f_ol,2))));

%% learning output weights for classic ring
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
result.baseline.z_ol = net.wout'*sim.r;

[x_cl,result.baseline.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));


figure;
plot(sim.psi,phase([1,1i]*result.baseline.z_ol)-sim.psi);


figure;plot(angle(result.baseline.z_cl)')


% [xZZ_cl,result.with_learningZZ.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',x,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
% 
% net.wfb=netbkp.wfb;
% [xFF_cl,result.with_learningFF.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',0.1*x,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
% 
% [xFF0_cl,result.with_learningFF0.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',ones(size(x)),...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
% 
[xFF2_cl,result.with_learningFF2.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',1e-4*ones(size(x)),...
'ol',0,...
'tmax',tmax,...
'dt',dt));
% 
% [xFF3_cl,result.with_learningFF3.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',x0a,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));