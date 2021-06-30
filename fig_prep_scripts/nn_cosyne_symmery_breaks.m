clear;
result=struct;
tmax=1000;
dt=1;

net.g=0;
net.N=1000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;
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

[x_cl,result.baseline.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));


%% Plot Wi-j
figure(1)
w=net.wfb*net.wout';
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);
theta=(0:net.N-1)./(net.N)*(2*pi)-pi;
plot(theta,wShift')
ylim([-0.2 0.2])
xlim([-pi pi])
box off
axis square
figure(2)
imagesc(w,[-0.02 0.02])


colorbar
box off
axis square


%% adding noise to 'classic' ring
result.with_noise={};
% g_vec=[0.01,1.0];
rng(1)
g_vec=[1.8];

trained_net={};
for gg=1:length(g_vec)
g=g_vec(gg);
net_with_g=net;
net = update_net_g(net,g);
x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
sim.r = net.phi(x_noi);
result.with_noise{gg}.z_ol = net.wout'*sim.r;
trained_net{gg}=net;
trained_net{gg}.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
[x_cl,result.with_noise{gg}.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));

trained_net{gg}.win=0;
[x_cl_trained,zzzzz]=fast_conv_to_fp_extended(trained_net{gg},[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));
end

%% Adding stimulus at 120 deg offset
net.win=net.wfb;
net0=net;
net0.W=0*net0.W;
trained_net{1}.win=trained_net{1}.wfb;
psi=40*pi/180;
epsilonmag=4e-2;

 
%  rot_result0=nn_simulate_rotation(sim,net0,x_cl(:,end),psi,epsilonmag);
 rot_result=nn_simulate_rotation(sim,net,x_cl(:,end),psi,epsilonmag);
 rot_result_trained=nn_simulate_rotation(sim,trained_net{1},x_cl(:,end),psi,epsilonmag);

%  nn_plot_rotation_results(sim,hp,rot_result0,psi,epsilonmag,'g=0 and no learning ');
 nn_plot_rotation_results(sim,hp,rot_result,psi,epsilonmag,...
     ['g=',num2str(trained_net{1}.g) ,' no learning ']);
 nn_plot_rotation_results(sim,hp,rot_result_trained,psi,epsilonmag,...
     ['g=',num2str(trained_net{1}.g) ,' with learning ']); 

%fix needed 
%  result.rotation=rot_result.rotation;
%  result.rotation_trained=rot_result.rotation;
% 
%  
% figure;plot(angle(result.rotation.zdelta(10:10:100,:))')
%% Plot Wi-j
w=net.wfb*net.wout'+net.W;
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);
theta=(0:net.N-1)./(net.N)*(2*pi)-pi;
figure(11)
plot(theta,wShift','.r')
hold all
plot(theta,mean(wShift,1),'k','LineWidth',3)
ylim([-0.2 0.2])
xlim([-pi pi])
box off
axis square
figure(22)
imagesc(w,[-0.04 0.04])
colorbar
box off
axis square
%%
% %% learning ring with g
% net = update_net_g(net,1.5);
% x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
% 
% sim.r = net.phi(x);
% result.with_learning={};
% M_vec=[5,20];
% for MM=1:length(M_vec)
% hp.M=M_vec(MM); 
% 
% pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
% regfac = hp.M*hp.alpha_reg*eye(length(pts));
% net.wout = (sim.r(:,pts)/...
%     (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
% 
% result.with_learning{MM}.z_ol = net.wout'*sim.r;
% [x_cl,result.with_learning{MM}.z_cl]=fast_conv_to_fp_extended(net,[],...
% struct('xinit',x,...
% 'ol',0,...
% 'tmax',tmax,...
% 'dt',dt));
% end

figure;
plot(sim.psi,phase([1,1i]*result.baseline.z_ol)-sim.psi);
hold on;
plot(sim.psi,phase([1,1i]*result.with_noise{1}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_noise{2}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi);
% plot(sim.psi,phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi);
legend('classic','small niose','large noise','M=5','M=20');

figure;plot(angle(result.baseline.z_cl(1:10:end,:))','.r')
ylim([-pi pi])
box off 
axis square
figure;plot(angle(result.with_noise{1}.z_cl(1:10:end,:))','.r')
ylim([-pi pi])
box off 
axis square
% figure;plot(angle(result.with_noise{2}.z_cl)')
% figure;plot(angle(result.with_learning{1}.z_cl)');
% figure;plot(angle(result.with_learning{2}.z_cl)');

figure;plot(angle(result.rotation.zdelta),'r');ylim([-pi pi])
box off 
axis square