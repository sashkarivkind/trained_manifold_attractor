clear;
saveFlag=0;
% saveFolder='/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/FiguresSubPlots/FigIntro/';
result=struct;
tmax=100;
dt=.1;

net.g=0;
net.N=4000;
net.phi=@(x)erf(x./sqrt(2));
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;
hp.M=40;
[hp, net, sim] = prep_network_param(hp, net,sim);

% genrate the initial conditions on the M points
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

% learning output weights for classic ring
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
% result.baseline.z_ol = net.wout'*sim.r;

% [result.baseline.x_cl,result.baseline.z_cl]=fast_conv_to_fp_extended(net,[],...
%     struct('xinit',x,...
%     'ol',0,...
%     'tmax',tmax,...
%     'dt',dt));
% adding noise to 'classic' ring
net0=net;
result.with_noise={};
% g_vec=[0.01,1.0];
rng(1)
g_vec=[0.01 0.2 1 1.5 1.8];
sigext=.1;
trained_net={};
for gg=1:length(g_vec)
    gg
    tic;
    g=g_vec(gg);
    net_with_g=net;
    net = update_net_g(net,g);
    x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
    sim.r = net.phi(x_noi);
%     result.with_noise{gg}.z_ol = net.wout'*sim.r;
    trained_net{gg}=net;
    trained_net{gg}.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
    trained_net{gg}.win=0;
    [x_cl_trained,result.trained{gg}.z_cl]=fast_conv_to_fp_extended_ou(trained_net{gg},[],...
        struct('xinit',x_noi,...
        'ol',0,...
        'tmax',tmax,...
        'dt',dt,'sigext',sigext));
    toc;
end
%% Have a look
for gg=1:length(g_vec)
    subplot(2,5,gg)

plot(dt:dt:tmax,ones(length(dt:dt:tmax),1).*angle(sim.f_ol(1,1:10:end)+1i.*sim.f_ol(2,1:10:end)),'--r');hold on; 
plot(dt:dt:tmax,(angle(result.trained{gg}.z_cl(1:10:end,:))),'.')
title(num2str(g_vec(gg)))
end

for gg=1:length(g_vec)
    subplot(2,5,gg+5)

plot(dt:dt:tmax,unwrap(angle(result.trained{gg}.z_cl(1:5:end,:))-angle(result.trained{gg}.z_cl(1:5:end,1))),'-')
title(num2str(g_vec(gg)))
ylim([-0.3 0.3])
end
%% Calculate diffusion coefficent
Dt_vec=dt:1*dt:tmax; % Dt=tau=10ms

ind=0;

for gg=1:length(g_vec)
angleZtest=unwrap(angle(result.trained{gg}.z_cl(1:end,:)));
    FFThet(gg,:,:)=diff(angleZtest(:,:),[],2);

for Dt=1:length(Dt_vec)
    %     ind=ind+1;
    % f(:,ind,:)=angle(Ztest(:,1+Dt,:))-angle(Ztest(:,1,:));
    % g(:,ind,:)=(angle(Ztest(:,1+Dt,:))-angle(Ztest(:,1,:))).^2;
    
    ind=ind+1;
%     fThet(gg,:,Dt)=angleZtest(:,Dt)-angleZtest(:,end/2);
%     gThet(gg,:,Dt)=(angleZtest(:,Dt)-angleZtest(:,end/2)).^2;

    fThet(gg,:,Dt)=angleZtest(:,Dt)-angleZtest(:,1);
%     fThet(gg,:,Dt)=diff(angleZtest(:,:),[],2);
    gThet(gg,:,Dt)=(angleZtest(:,Dt)-angleZtest(:,1)).^2;
    
    
%     fMod(ind,:,:)=modZtest(1+Dt,:,:)-modZtest(1,:,:);
%     gMod(ind,:,:)=(modZtest(1+Dt,:,:)-modZtest(1,:,:)).^2;
    
end
end
figure
subplot(2,1,1)
plot(Dt_vec,squeeze(mean(gThet,2)))
legend('g=0.01','g=0.2','g=1','g=1.5','g=1.8')
xlabel('time (tau)'); ylabel('variance: [(\psi(end/2+dt)-\psi(end/2))^2]_{angels}/dt')
title('Diffusion coeff: Var of psi Vs time- linear in diffusion process')

subplot(2,1,2)
plot(Dt_vec,squeeze(mean(fThet,2))./dt)
legend('g=0.01','g=0.2','g=1','g=1.5','g=1.8')
xlabel('time (tau)'); ylabel('mean: [(\psi(end/2+dt)-\psi(end/2))]_{angels}/dt')
title('Drift: mean of psi Vs time- constant not on a ring')

figure
plot(angleZtest(:,5),squeeze(mean(FFThet(3,:,:),3))./dt,'.')