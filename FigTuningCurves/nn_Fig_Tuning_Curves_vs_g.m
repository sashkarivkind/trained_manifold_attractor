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

sim=struct;hp.A=1.5;
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
% g_vec=[0.01 0.1:0.1:1.4];
g_vec=[0.01 0.1:0.1:1.8];

trained_net={};
for gg=1:length(g_vec)
    g=g_vec(gg)
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

%%
% % Plot Tuning curves
% % Classic
% xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
% figMemTuning=figure(23);plot(xaxis_deg-pi,(net.phi(result.baseline.x_cl(1:200:1000,:)))');
% xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])
% % Learned
% xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
% figMemTuning2=figure(24);plot(xaxis_deg-pi,(net.phi(result.with_learning{1}.x_cl(1:100:1000,:)))');
% xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])

%% Calculate statistics
% close all
Dpr=[];
fuBase=fft(net.phi(result.baseline.x_cl(1:end,:)),[],2);

for gg=1:length(g_vec)
fuRand=fft(net.phi(result.with_learning{gg}.x_cl(1:end,:)),[],2);
f1_stat(gg,:)=abs(fuRand(:,2))./size(fuRand,2);
mf1_stat(gg)=mean(f1_stat(gg,:));
sdf1_stat(gg)=std(f1_stat(gg,:));


r_o1=(net.phi(result.with_learning{gg}.x_cl(1:end,:))+1)./2;
fuRand_01=fft(r_o1,[],2);
OSI_stat_01(gg,:)=1-abs(fuRand_01(:,2))./abs(fuRand_01(:,1));
mOSI_stat(gg)=mean(OSI_stat_01(gg,:));
sdOSI_stat(gg)=std(OSI_stat_01(gg,:));

figure(51)
subplot (5,4,gg)
hist(f1_stat(gg,:),50)
xlim([0 1])
title(num2str(g_vec(gg)))
xlabel('f_1'); ylabel('counts')
box off

figure(52)
subplot (5,4,gg)
hist(OSI_stat_01(gg,:),50)
xlim([0 1])
title(num2str(g_vec(gg)))
xlabel('OSI'); ylabel('counts')
box off

% For PR
r=net.phi(result.with_learning{gg}.x_cl(1:end,:));
C=r'*r; lambda=eig(C,'Vector');
lambdaSort=sort(lambda,'descend');

% lambdasum2= sum(lambdaSort(1:2*hp.M)).^2;
% lambda2sum= sum( (lambdaSort(1:2*hp.M)).^2  );
lambdasum2= sum(lambdaSort(1:6)).^2;
lambda2sum= sum( (lambdaSort(1:6)).^2  );
% lambdasum2= sum(lambda).^2;
% lambda2sum= sum( (lambda).^2  );
Dpr(gg)=lambdasum2./lambda2sum;
end
figure(100)
plot(g_vec,sdf1_stat,'o');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Tuning curves diversity (SD (f_1) )');

figure(101)
plot(g_vec,sdOSI_stat,'o');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Tuning curves diversity (SD (OSI) )');


figure(200)
plot(g_vec,Dpr,'o');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Participation ratio ((\Sigma\lambda)^2/(\Sigma \lambda^2) )');

figure(300)
plot(g_vec,mf1_stat,'o');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Tuning curves mean (mean (f_1) )');

figure(301)
plot(g_vec,mOSI_stat,'o');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Tuning curves mean (mean (OSI) )');

%% Generate MF estimate of te sd(f^1_i)
% N=100000;

N=1000;
tpN=[0:(N-1)]/N*2*pi;tpN=tpN';
Amp=hp.A; res=160;
psi_vec=[0:(res-1)]/res*2*pi;
for gg=1:length(g_vec)
    gg
    g=g_vec(gg);
    % Solve mean field to get correaltion function
    [mu, sigxa,cxa,ca]=...
        solveFixedPointRing_Erf(net.phi,0,g,Amp,1,1,100,psi_vec);
    sigx=sqrt(sigxa.^2-Amp^2/2); cx=cxa-Amp.^2.*ca; Cr=cx./g^2;
    cfu=fft(Cr);
    cfu=2*real(cfu)/length(cfu); %cosine fourier coefficients
    hmax=10;
    c=g^2*cfu(2:1:hmax);
    
    % Get MF result from sd(a_1)/g ????
    reCo=real(sqrt(real(c(1:2:end)))); %real odd fourier coefficients
    sdf1_stat_MF_test(gg)=reCo(1)/g;
    
    % Get MF result from numerically solving the integral
    ai=randn(N,length(reCo))*diag(reCo); %cos coefficient
    bi=randn(N,length(reCo))*diag(reCo); %sin coefficient
    for ii=1:length(psi_vec)
        harmonics=psi_vec(ii)'*[1:2:length(c)];
        x1_approx=sum((ai.*cos(harmonics))...
            +(bi.*sin(harmonics)),2);
        result.with_learning{gg}.MF_TC(:,ii)=net.phi(Amp*cos(tpN-psi_vec(ii))+x1_approx);
    end
    
    % Calculate f_1
    fuRandMF=fft(result.with_learning{gg}.MF_TC,[],2);
    f1_stat_MF(gg,:)=abs(fuRandMF(:,2))./size(fuRandMF,2);
    sdf1_stat_MF(gg)=std(f1_stat_MF(gg,:));
    % Plot the statistics
    figure(33)
    subplot (4,4,gg)
    hist(f1_stat_MF(gg,:),50)
    xlim([0 1])
    title(num2str(g_vec(gg)))
    xlabel('f_1'); ylabel('counts')
    box off
    
    
    % Calculate CirVar
    r_o1_MF=(result.with_learning{gg}.MF_TC+1)./2;
    fuRand_01_MF=fft(r_o1_MF,[],2);
    OSI_stat_01_MF(gg,:)=1-abs(fuRand_01_MF(:,2))./abs(fuRand_01_MF(:,1));
    mOSI_stat_MF(gg)=mean(OSI_stat_01_MF(gg,:));
    sdOSI_stat_MF(gg)=std(OSI_stat_01_MF(gg,:));
    % Plot the statistics
    figure(34)
    subplot (4,4,gg)
    hist(OSI_stat_01_MF(gg,:),50)
    xlim([0 1])
    title(num2str(g_vec(gg)))
    xlabel('CircVar'); ylabel('counts')
    box off
    
end

figure(100)
plot(g_vec,sdf1_stat_MF,'-k');hold on
% plot(g_vec,sdf1_stat_MF_test,'-r');hold on
%% Plot Tuning curves example with high and low SI
gg=15 %0.2;
% id=[20 21 22];
id=[520 182 521 609];
% id=1:1000;
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning2=figure(24);
% plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(end/2+id,:)))');
plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(id,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])


fuRand=fft(net.phi(result.with_learning{gg}.x_cl(id,:)),[],2);

f1_stat_id=abs(fuRand(:,2))./size(fuRand,2)
    title(['g=1.4 f_1=',num2str(f1_stat_id')])

[m mm]=sort(f1_stat_id);
%% Calculate participation ratio

