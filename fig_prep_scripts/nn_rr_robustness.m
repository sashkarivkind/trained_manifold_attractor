clear;
result=struct;
tmax=1000;
dt=0.1;
input_mag=0.01;
input_noi=3e-3;
omega=0; %% 0 for white noise
net.g=0;
net.N=4000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;

rng(1)
g_vec=[0.01,0.5:0.5:1.5, 1.8];
GrrCL0 =[    2.3215
    1.7491
    0.8902
    0.3824
    0.2674];
trained_net={};
for gg=1:length(g_vec)
    tic;
g=g_vec(gg)
[hp, net, sim] = prep_network_param(hp, net,sim);
net = update_net_g(net,g);
xol = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
sim.r = net.phi(xol);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
trained_net{gg}=net;
trained_net{gg}.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
[x_cl,result.with_noise{gg}.z_cl]=fast_conv_to_fp_extended(trained_net{gg},[],...
struct('xinit',xol,...
'ol',0,...
'tmax',tmax,...
'dt',dt));

% trained_net{gg}.win=0;
% [x_cl_trained,zzzzz]=fast_conv_to_fp_extended(trained_net{gg},[],...
% struct('xinit',x,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
trained_net{gg}.win=trained_net{gg}.wfb;
rrnoi_result_trained{gg}=nn_simulate_noisy_stimulus(sim,trained_net{gg},x_cl(:,1),input_mag,input_noi,omega);
toc;
end
%%
figure; 
for gg=1:length(g_vec)
plot(abs(rrnoi_result_trained{gg}.noisy.zdelta));
hold on;
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
plot(hp.A+GrrCL0(gg)*input_mag*ones(size(rrnoi_result_trained{gg}.noisy.zdelta)),'--');
% plot(angle(rrnoi_result_trained{gg}.noisy.zdelta));
% zabs_bar(gg)=mean(abs(rrnoi_result_trained{gg}.noisy.zdelta(end-2000:end)));
% zabs_std(gg)=std(abs(rrnoi_result_trained{gg}.noisy.zdelta(end-2000:end)));
hold on;
end
legend('0.01','0.5','1','1.5','1.8')