clear;
saveFlag=0;
saveFolder='/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/FiguresSubPlots/FigLearnExample/';
result=struct;
tmax=1000;
dt=1;
rng(2)
net.g=0;
net.N=1000;
net.phi=@(x)erf(x./sqrt(2));
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;
hp.M=40;
hp.A=1.2;
[hp, net, sim] = prep_network_param(hp, net,sim);
hp
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

% g_vec=[0 0.5 1 1.2 1.5 1.7];
g_vec=[0.5 1.2 1.7];
theta=(0:net.N-1)./(net.N)*(2*pi)-pi;


trained_net={};
for gg=1:length(g_vec)

    g=g_vec(gg);
    net_with_g=net;
    net = update_net_g(net,g);
    x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
    sim.r = net.phi(x_noi);
%     result.with_noise{gg}.z_ol = net.wout'*sim.r;
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






w=trained_net{gg}.wfb*trained_net{gg}.wout'+trained_net{gg}.W;
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);
figSym3=figure(1000*(gg-1)+13);plot(theta,wShift(1:50:end,:)','.r');hold all;plot(theta,mean(wShift,1),'k','LineWidth',3)
ylim([-0.1 0.1]); xlim([-pi pi]);box off; axis square
% figMat3=figure(1000*(gg-1)+14);imagesc(w,[-0.1 0.1]);colorbar; box off; axis square

w=trained_net{gg}.wfb*trained_net{gg}.wout';
wShift=w;
for ii=1:net.N
    wShift(ii,:)=circshift(w(ii,:),-(ii-1));
end
wShift=fftshift(wShift);
figSym3=figure(1000*(gg-1)+113);plot(theta,wShift(1:50:end,:)','.r');hold all;plot(theta,mean(wShift,1),'k','LineWidth',3)
ylim([-0.01 0.01]); xlim([-pi pi]);box off; axis square
% figMat3=figure(1000*(gg-1)+114);imagesc(w,[-0.01 0.01]);colorbar; box off; axis square

% Learned
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning2=figure(1000*(gg-1)+24);plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(1:100:1000,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])


end
%%
%
% figMemClassic=figure(21);plot(angle(result.baseline.z_cl(1:10:end,:))','.r')
% ylim([-pi pi]);box off ;axis square
% 
% figMemNoise=figure(22);plot(angle(result.with_noise{gg}.z_cl(1:10:end,:))','.r');
% ylim([-pi pi])
% box off ;axis square
% 
% figMemLearn=figure(23);plot(angle(result.with_learning{gg}.z_cl(1:10:end,:))','.r');
% ylim([-pi pi])
% box off ;axis square

% Save Tuning curves

