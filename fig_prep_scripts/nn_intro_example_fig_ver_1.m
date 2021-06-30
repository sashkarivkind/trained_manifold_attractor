clear;
saveFlag=0;
saveFolder='/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/FiguresSubPlots/FigLearnExample/';
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
g_vec=[1.2];

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
    
    
    trained_net{gg}.win=0;%%%%%%%%%%%%%%%%%%%%%%%%%% result.training....
    result.with_learning{gg}.z_ol = net.wout'*sim.r;

    [result.with_learning{gg}.x_cl,result.with_learning{gg}.z_cl]=fast_conv_to_fp_extended(trained_net{gg},[],...
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

psi_traget_vec=[-2*pi/3 0 pi/2];
%
for psi=psi_traget_vec
    
    
    epsilonmag=4e-2;
    % Classic ring
    rot_result0=nn_simulate_rotation(sim,net0,x_cl(:,end),psi,epsilonmag);
    figRotDecdoerErrorClassic=figure(100);
    xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
    res=length(xaxis_deg);
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result0.rotation.zdelta_ol)-sim.psi,res/2  )  )
    hold on
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result0.rotation.zdelta_ol_ref)-sim.psi,res/2  )  )
    xlim([-pi pi]); box off
    ylim([-0.1 0.1]);set(gca,'Xtick',[-pi:pi/4:pi])
    
    figRotTrajClassic=figure(102);
    plot(angle(rot_result0.rotation.zdelta))
    hold on
    plot(psi.*ones(length(rot_result0.rotation.zdelta),1),'k')
    ylim([-pi pi]); box off ; set(gca,'Ytick',[-pi:pi/4:pi])
    
    % Ring+noise
    rot_result=nn_simulate_rotation(sim,net,x_cl(:,end),psi,epsilonmag);
    figRotDecdoerErrorNoise=figure(101);
    xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result.rotation.zdelta_ol)-sim.psi,res/2  )  )
    hold on
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result.rotation.zdelta_ol_ref)-sim.psi,res/2  )  )
    plot(xaxis_deg-pi,zeros(length(rot_result.rotation.zdelta_ol_ref),1),'k')
    xlim([-pi pi]); box off
    ylim([-0.1 0.1]);set(gca,'Xtick',[-pi:pi/4:pi])
    
    figRotTrajNoise=figure(103);
    plot(angle(rot_result.rotation.zdelta))
    hold on
    plot(angle(rot_result.rotation.zdelta))
    plot(psi.*ones(length(rot_result.rotation.zdelta),1),'k')
    ylim([-pi pi]); box off ; set(gca,'Ytick',[-pi:pi/4:pi])
    
    % Learned Ring+noise
    rot_result_trained=nn_simulate_rotation(sim,trained_net{1},x_cl(:,end),psi,epsilonmag);
    figRotDecdoerErrorLearn=figure(104);
    xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result_trained.rotation.zdelta_ol)-sim.psi,res/2  )  )
    hold on
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result_trained.rotation.zdelta_ol_ref)-sim.psi,res/2  )  )
    plot(xaxis_deg-pi,zeros(length(rot_result_trained.rotation.zdelta_ol_ref),1),'k')
    xlim([-pi pi]); box off
    ylim([-0.1 0.1]);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    figRotDecdoerErrorLearnZoom=figure(1004);
    xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
    plot(xaxis_deg-pi,circshift(phase([1,1i]*rot_result_trained.rotation.zdelta_ol_ref)-sim.psi,res/2  )  )
    hold on
    plot(xaxis_deg-pi,zeros(length(rot_result_trained.rotation.zdelta_ol_ref),1),'k')
    xlim([-pi pi]); box off
    ylim([-0.00001 0.00001]);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    figRotTrajLearn=figure(105);
    plot(angle(rot_result_trained.rotation.zdelta))
    hold on
    plot(angle(rot_result_trained.rotation.zdelta))
    plot(psi.*ones(length(rot_result_trained.rotation.zdelta),1),'k')
    ylim([-pi pi]); box off ; set(gca,'Ytick',[-pi:pi/4:pi])
    
end

if saveFlag
    saveas(figRotDecdoerErrorClassic,[saveFolder 'RotDecoderErClassic.pdf'])
    saveas(figRotDecdoerErrorNoise,[saveFolder 'RotDecoderErNoise.pdf'])
    saveas(figRotTrajClassic,[saveFolder 'figRotTrajClassic.pdf'])
    saveas(figRotTrajNoise,[saveFolder 'figRotTrajNoise.pdf'])
    saveas(figRotTrajLearn,[saveFolder 'figRotTrajLearn.pdf'])
    saveas(figRotDecdoerErrorLearn,[saveFolder 'RotDecoderErLearn.pdf'])
    saveas(figRotDecdoerErrorLearnZoom,[saveFolder 'figRotDecdoerErrorLearnZoom.pdf'])
end
result.rotation0=rot_result0.rotation;
result.rotation=rot_result.rotation;
%  result.rotation_trained=rot_result.rotation;

% Plot W matrix
%
w=net.wfb*net.wout';
wShift=w;
theta=(0:net.N-1)./(net.N)*(2*pi)-pi;

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
figMat2=figure(12);imagesc(w,[-0.05 0.05]);colorbar; box off; axis square
%
w=trained_net{1}.wfb*trained_net{1}.wout'+trained_net{1}.W;
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);

figSym3=figure(13);plot(theta,wShift(1:50:end,:)','.r');hold all;plot(theta,mean(wShift,1),'k','LineWidth',3)
ylim([-0.1 0.1]); xlim([-pi pi]);box off; axis square
figMat3=figure(14);imagesc(w,[-0.1 0.1]);colorbar; box off; axis square


if saveFlag
    saveas(figSym1,[saveFolder 'ClassicSym.pdf'])
    saveas(figMat1,[saveFolder 'ClassicMat.pdf'])
    saveas(figSym2,[saveFolder 'NoiseSym.pdf'])
    saveas(figMat2,[saveFolder 'NoiseMat.pdf'])
     saveas(figSym3,[saveFolder 'LearnSym.pdf'])
    saveas(figMat3,[saveFolder 'LearnMat.pdf'])
end
%
figMemClassic=figure(21);plot(angle(result.baseline.z_cl(1:10:end,:))','.r')
ylim([-pi pi]);box off ;axis square

figMemNoise=figure(22);plot(angle(result.with_noise{1}.z_cl(1:10:end,:))','.r');
ylim([-pi pi])
box off ;axis square

figMemLearn=figure(23);plot(angle(result.with_learning{1}.z_cl(1:10:end,:))','.r');
ylim([-pi pi])
box off ;axis square

if saveFlag
    saveas(figMemClassic,[saveFolder 'MemoryClassic.pdf'])
    saveas(figMemNoise,[saveFolder 'MemoryNoise.pdf'])
    saveas(figMemLearn,[saveFolder 'figMemLearn.pdf'])
end
% Save Tuning curves
% Classic
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning=figure(23);plot(xaxis_deg-pi,(net.phi(result.baseline.x_cl(1:200:1000,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])
% Learned
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning2=figure(24);plot(xaxis_deg-pi,(net.phi(result.with_learning{1}.x_cl(1:100:1000,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])


if saveFlag
    saveas(figMemTuning,[saveFolder 'TunincurvesClassic.pdf'])
    saveas(figMemTuning2,[saveFolder 'TunincurvesLearn.pdf'])
end