clear;
saveFlag=0;
saveFolder='/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/FiguresSubPlots/FigIntro/';
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

% learning output weights for classic ring
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
result.baseline.z_ol = net.wout'*sim.r;

[result.baseline.x_cl,result.baseline.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));
% adding noise to 'classic' ring
result.with_noise={};
% g_vec=[0.01,1.0];
rng(1)
g_vec=[.5];

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

% Adding stimulus at 120 deg offset
net.win=net.wfb;
net0=net;
net0.W=0*net0.W;
trained_net{1}.win=trained_net{1}.wfb;
% psi=40*pi/180;

psi=2*pi/4;
epsilonmag=4e-2;

rot_result0=nn_simulate_rotation(sim,net0,x_cl(:,end),psi,epsilonmag);
figRotDecdoerErrorClassic=figure(100);
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
plot(xaxis_deg,phase([1,1i]*rot_result0.rotation.zdelta_ol)-sim.psi)
hold on
plot(xaxis_deg,phase([1,1i]*rot_result0.rotation.zdelta_ol_ref)-sim.psi)
xlim([0,2*pi]); box off
ylim([-0.06 0.06]);set(gca,'Xtick',[0:pi/4:2*pi])

rot_result=nn_simulate_rotation(sim,net,x_cl(:,end),psi,epsilonmag);
figRotDecdoerErrorNoise=figure(101);
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
plot(xaxis_deg,phase([1,1i]*rot_result.rotation.zdelta_ol)-sim.psi)
hold on
plot(xaxis_deg,phase([1,1i]*rot_result.rotation.zdelta_ol_ref)-sim.psi);
plot(xaxis_deg,zeros(length(rot_result.rotation.zdelta_ol_ref),1),'k')
xlim([0,2*pi]); box off
ylim([-0.06 0.06]);set(gca,'Xtick',[0:pi/4:2*pi])

figRotTraj=figure(102);
plot(angle(rot_result0.rotation.zdelta))
hold on
plot(angle(rot_result.rotation.zdelta))
plot(psi.*ones(length(rot_result.rotation.zdelta),1),'k')
ylim([-pi pi]); box off ; set(gca,'Ytick',[-pi:pi/4:pi])

if saveFlag
        saveas(figRotDecdoerErrorClassic,[saveFolder 'RotDecoderErClassic.pdf'])
        saveas(figRotDecdoerErrorNoise,[saveFolder 'RotDecoderErNoise.pdf'])
        saveas(figRotTraj,[saveFolder 'RotTraj.pdf'])
end
result.rotation0=rot_result0.rotation;
result.rotation=rot_result.rotation;
%  result.rotation_trained=rot_result.rotation;

% Plot W matrix
%%
w=net.wfb*net.wout';
wShift=w;

figSym1=figure(1);plot(theta,fftshift(w(:,1)));ylim([-0.05 0.05]); xlim([-pi pi]);box off; axis square
figMat1=figure(2);imagesc(w,[-0.02 0.02]);colorbar; box off; axis square

w=net.wfb*net.wout'+net.W;
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);

figSym2=figure(11);plot(theta,wShift(1:50:end,:)','.r');hold all;plot(theta,mean(wShift,1),'k','LineWidth',3)
ylim([-0.05 0.05]); xlim([-pi pi]);box off; axis square
figMat2=figure(12);imagesc(w,[-0.04 0.04]);colorbar; box off; axis square

if saveFlag
        saveas(figSym1,[saveFolder 'ClassicSym.pdf'])
        saveas(figMat1,[saveFolder 'ClassicMat.pdf'])
        saveas(figSym2,[saveFolder 'NoiseSym.pdf'])
        saveas(figMat2,[saveFolder 'NoiseMat.pdf'])
end
%%
figMemClassic=figure(21);plot(angle(result.baseline.z_cl(1:10:end,:))','.r')
ylim([-pi pi]);box off ;axis square

figMemNoise=figure(22);plot(angle(result.with_noise{1}.z_cl(1:10:end,:))','.r');
ylim([-pi pi])
box off ;axis square

if saveFlag
        saveas(figMemClassic,[saveFolder 'MemoryClassic.pdf'])
        saveas(figMemNoise,[saveFolder 'MemoryNoise.pdf'])
end
% Save Tuning curves
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning=figure(23);plot(xaxis_deg-pi,(net.phi(result.baseline.x_cl(1:200:1000,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])
if saveFlag
        saveas(figMemTuning,[saveFolder 'TunincurvesClassic.pdf'])
end