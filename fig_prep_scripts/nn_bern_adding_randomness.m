clear;
result=struct;
tmax=1000;
dt=1;

net.g=0;
net.N=1000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);

[hp, net, sim] = prep_network_param(struct, net,struct);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%% learning output weights for classic ring
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
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


%% adding noise to 'classic' ring
result.with_noise={};
g_vec=[0.01,1.0];
for gg=1:length(g_vec)
g=g_vec(gg);
net_with_g=net;
net = update_net_g(net,g);
x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
sim.r = net.phi(x_noi);
result.with_noise{gg}.z_ol = net.wout'*sim.r;

[x_cl,result.with_noise{gg}.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));
end

%% learning ring with g
net = update_net_g(net,1.5);
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

sim.r = net.phi(x);
result.with_learning={};
M_vec=[5,20];
for MM=1:length(M_vec)
hp.M=M_vec(MM);

pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

result.with_learning{MM}.z_ol = net.wout'*sim.r;
[x_cl,result.with_learning{MM}.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));
end

figure;
plot(sim.psi,phase([1,1i]*result.baseline.z_ol)-sim.psi);
hold on;
plot(sim.psi,phase([1,1i]*result.with_noise{1}.z_ol)-sim.psi);
plot(sim.psi,phase([1,1i]*result.with_noise{2}.z_ol)-sim.psi);
plot(sim.psi,phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi);
plot(sim.psi,phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi);
legend('classic','small niose','large noise','M=5','M=20');

figure;plot(angle(result.baseline.z_cl)')
figure;plot(angle(result.with_noise{1}.z_cl)')
figure;plot(angle(result.with_noise{2}.z_cl)')
figure;plot(angle(result.with_learning{1}.z_cl)');
figure;plot(angle(result.with_learning{2}.z_cl)');

