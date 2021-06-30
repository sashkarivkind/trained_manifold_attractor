% clear;
result=struct;
tmax=1000;
%%measuring drift speed
dt=0.1;
tmax_drift=3;
velocity_sample_timestep=7:8;
%%net settings
net.g=0;
net.N=2000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;
hp.A=1.2;
[hp, net, sim] = prep_network_param(hp, net,sim);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%% learning output weights for classic ring
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%% obtaining open loop output
result.baseline.z_ol = net.wout'*sim.r;

% [x_cl,result.baseline.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',x,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
% 



%% adding noise to 'classic' ring
result.with_noise={};
% g_vec=[0.01,1.0];
rng(87)
g_vec=[0.5];

trained_net={};
for gg=1:length(g_vec)
g=g_vec(gg);
net_with_g=net;
net = update_net_g(net,g);
x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
sim.r = net.phi(x_noi);
result.with_noise{gg}.z_ol = net.wout'*sim.r;
end
%%
delta_psi=phase([1,1i]*result.with_noise{1}.z_ol)-sim.psi;
delta_rho=abs([1,1i]*result.with_noise{1}.z_ol)-hp.A;
rmse=sqrt(mean(delta_psi.^2))
eDC=mean(delta_psi)
rmseNoDc=sqrt(mean(delta_psi.^2)-eDC^2)
%%
[x_cl,result.with_noise{gg}.z_cl]=fast_conv_to_fp_extended(net,[],...
    struct('xinit',x,...
    'ol',0,...
    'tmax',tmax_drift,...
    'dt',dt,'save_neurons',1));
%%
        v=1/dt*mean(diff(mod(angle(result.with_noise{gg}.z_cl(:,velocity_sample_timestep)),2*pi),[],2),2);
        figure; plot(v',(phase([1,1i]*result.with_noise{1}.z_ol)-sim.psi),'x')
        
%%

figure; 
plot(delta_psi,v,'x')
hold on;
plot([-1,1],[-1,1],'k--');
xlim(1.5*[-1,1]*max(abs(delta_psi)));
ylim(1.5*[-1,1]*max(abs(delta_psi)));
xlabel('\Delta');
ylabel('v_{drift}')
%%
figure;
plot(sim.psi,phase([1,1i]*result.baseline.z_ol)-sim.psi);
hold on;
plot(sim.psi,phase([1,1i]*result.with_noise{1}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_noise{2}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi);
legend('classic','small niose','large noise','M=5','M=20');


% figure;plot(angle(result.baseline.z_cl(1:10:end,:))','.r')
% ylim([-pi pi])
% box off 
% axis square
% figure;plot(angle(result.with_noise{1}.z_cl(1:10:end,:))','.r')
% ylim([-pi pi])
% box off 
% axis square

%%
% figure; 
% subplot(2,1,1);
% plot(angle(result.with_noise{1}.z_cl).');
% subplot(2,1,2);
% plot(abs(result.with_noise{1}.z_cl).');
