% clear;
% close all
tau=1;
% Amp=0.033;N=3;J_0=0;
Amp=0.2;
TrainParam.UniformFlag=0;

% Kappa controls the width of the bump.
kappa=2/Amp;
NetParm.a=6;NetParm.b=1;NetParm.h=1;
M=50;
mm=1;
J_0=0.1;
% 
% NetParm.JType='ER';
% phi=@(x)x.*(x>0);
% phi_prim=@(x) 1.0*(x>0);
% N=2000;t_end=500; K=400; J_ext=1.*J_0;del_t=tau/20;

% NetParm.JType='ER';
% phi=@(x)x.^2.5.*(x>0);
% phi_prim=@(x) 2.5.*x.^(1.5).*(x>0);
% N=2000;t_end=500; K=400; J_ext=0.2.*J_0;del_t=tau/20;

    NetParm.JType='Gaussian';
    phi=@(x)tanh(x);    
    phi_prim=@(x)(1-tanh(x).^2);
% phi=@(x)erf(x./sqrt(2));
% phi_prim =  @(x) exp(-x.^2/2)/sqrt(pi/2);

N=1000;t_end=500; K=100; J_ext=0;del_t=tau;
% This is the new non linear-function. If you want to go back to
% cos(theta-psi) you should take the identity function
% G=@(A) 3.*(exp(kappa.*A )./(2*pi.*besseli(0,kappa))-1/(2*pi)); %
G=@(A) A; %identity
NetParm.G=G;
%Regularization
TrainParam.alpha=0e-10;

% TrainParam.alpha=5e-4;



seed=31;
seedUnlearned=20;

ii=0;
gUnlearned=0.0;
for S=[0 ]
    NetParm.S=S;
    %     NetParm.J_extCL=J_extCL;
    %     J_0=1;
    ii=ii+1;
    % Parameters for an inhibitory balanced network
    % For ReLu the transition is at J_0=sqrt(2), independent of J_ext
    % For other non-linearites- see Harish and Hansel 2015 or (Kadmon 2015)
    % mean rate in the network is given by the balanced equation:rate=J_ext/J_0
    
    
    %
    
    
    
    NetParm.J_extCL=J_ext;
    
    AmpC=1;
    AmpS=1; % to make an ellipse change this parameter
    
    NetParmUnlearned.JType=NetParm.JType;
    NetParmUnlearned.J_0=gUnlearned;
    NetParmUnlearned.N=N; NetParmUnlearned.K=K;
    NetParmUnlearned.seed=seedUnlearned;
    
    % Least square or PopVec decoder
    NetParm.phi=phi; NetParm.J_0=J_0;NetParm.J_ext=J_ext;
    NetParm.N=N;NetParm.K=K; NetParm.seed=seed; NetParm.tau=tau;
    
    TrainParam.AmpC=AmpC;TrainParam.AmpS=AmpS;
    TrainParam.Amp=Amp; TrainParam.t_end=400; TestParam.t_end=t_end;
    % Flags. Should all be at 0!
    PlotFlagTrain=0;
    PlotFlagTest=0;
    LearnPopVecFlag=0;
    % For learning with Chaos, take the mean rate instead of the final point
    ChaosFlag=0;
    TestParam.Mtest=M;
    TrainParam.Mtrain=M;
    M=TrainParam.Mtrain;
    
    TrainParam.del_t=del_t;
    TestParam.del_t=del_t;
    
    TestParam.SigTest=0.0;
    TestParam.sigIC=0.0;
    TestParam.Trials=1;
    
    TestParam.seedIC=1;
    tic;
    [rmsTrain, RmsTest, Ztest,Ysave,wfffb,normWout]=RunRingNoisySelfsyn(NetParm,NetParmUnlearned,TrainParam,TestParam,PlotFlagTrain,PlotFlagTest,ChaosFlag,LearnPopVecFlag);
    toc;
    
    nWout(mm)=normWout;
    
    figure(390)
    subplot 311
    plot(rmsTrain.totErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Error in Amp+Angle')
    subplot 312
    plot(rmsTrain.angErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Error in Angle')
    subplot 313
    plot(rmsTrain.ampErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Error in Amplitude')
    
    figure(391)
    subplot 311
    plot(RmsTest.totErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Test Error in Amp+Angle')
    subplot 312
    plot(RmsTest.angErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Test Error in Angle')
    subplot 313
    plot(RmsTest.ampErrTrain)%,'Color',col(mm,:));
    hold on
    xlabel('time')
    ylabel('Test Error in Amplitude')
    %
    figure(J_0*1000+ii*100+seed)
    % dt=tau/20;
    % tCut=100*round(tau/dt);
    tCut=1;
    ZtestCut=Ztest;
    ZtestCut(1:tCut,:,:)=[];
    angleZtest=unwrap(angle(ZtestCut),[],2);
    modZtest=abs(ZtestCut);
    subplot(2,1,1)
    plot(angleZtest)
    subplot(2,1,2)
    plot(modZtest)
    
    
    figure(300+ii)
    plot(real(Ztest(end,:)),imag(Ztest(end,:)),'o')
    hold on
    
    a=NetParm.a;b=NetParm.b;h=NetParm.h;
    psi_vec=2*pi*(0:TrainParam.Mtrain*2-1)/(TrainParam.Mtrain*2);
    f=Amp*[(a-b)*cos(psi_vec)+h*cos((a-b)/b*psi_vec); (a-b)*sin(psi_vec)-h*sin((a-b)/b*psi_vec)];
    if TrainParam.UniformFlag
        f = curvspace(f',length(psi_vec));
        f=f';
    end
    plot(f(1,:),f(2,:),'xr')
end

J=GenerateJ(NetParm);
%
for ii=1:M
    ii
    z=phi_prim(Ysave(:,ii,1));zsort=sort(eig(diag(z)*(J+wfffb)));
    if (strcmp(NetParm.JType,'ER'))
        gr1(ii)=zsort(end-1);gr2(ii)=zsort(end-2);
        
    else
        gr1(ii)=zsort(end);gr2(ii)=zsort(end-1);
        
    end
end
figure(12)
plot(abs(gr1),'r')
hold all
plot(abs(gr2),'b')
legend('ev1','ev2')

ii=3
z=phi_prim(Ysave(:,ii,1));zsort=sort(eig(-eye(N)+diag(z)*(J+wfffb)));
figure;plot(zsort,'o')
%%
pt=1;
ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
ee_cl1=eig((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt))));
ee_cl2=eig((net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(x(:,pt))));
ee_clg0=eig((0*net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));

figure;
plot(ee_cl,'*');
hold on;
plot(ee_clg0,'d');
plot(ee_cl1,'o');
plot(ee_cl2,'s');
%%
for ii=1:2:M
z=phi_prim(Ysave(:,ii,1));zsort=sort(eig(-eye(N)+diag(z)*(J+wfffb)));
figure;plot(zsort,'o')
end
%%
figure(33)

for tt=1:10:(size(Ztest,1)-1)

    plot(f(1,:),f(2,:),'-r')
    hold all
    box off
    axis square
    plot(real(Ztest(tt,:)),imag(Ztest(tt,:)),'o')
%     plot(angleZtest(tt),'o')   
    xlim([-4 4])
    ylim([-4 4])
    title(num2str(tt))
    pause(.1)
    hold off
end
%%
figure(33)
vidfile = VideoWriter('HypoM50g1_0_h_0_5_Amp_3_alpha_1e-10.mp4','MPEG-4');
open(vidfile);
for tt=1:10:(size(Ztest,1)-1)
    subplot(2,1,1)
    plot(abs(gr1),'r')
    hold all
    plot(abs(gr2),'b')
    xlim([0 50])
    axis square
    subplot(2,1,2)
    plot(f(1,:),f(2,:),'-r')
    set(gcf,'position',[742   508   378   440])

    hold all
    box off
    axis square
    plot(real(Ztest(tt,:)),imag(Ztest(tt,:)),'o')
%     plot(angleZtest(tt),'o')   
    xlim([-4 4])
    ylim([-4 4])
    title(num2str(tt))
%     pause(.1)
    hold off
        F(ii) = getframe(gcf);
    writeVideo(vidfile,F(ii));
end
close(vidfile)
%% Save video
vidfile = VideoWriter('BumpN50M6Inp0_1_g0_7.mp4','MPEG-4');
open(vidfile);

for ii=1:size(TuningCurves,2)
    ii
    plot(TuningCurves(:,ii),'LineWidth',3)
    set(gcf,'position',[806   731   314   217])
    box off
    ylim([-1 1])
    F(ii) = getframe(gcf);
    writeVideo(vidfile,F(ii));
end
close(vidfile)
%  box off
%%
Plot points on the attractor during the test
for ttt=1:10:size(Ysave,3)
    
    figure(300+ii)
    plot(real(Ztest(ttt,:)),imag(Ztest(ttt,:)),'o')
    hold on
    
    a=NetParm.a;b=NetParm.b;h=NetParm.h;
    psi_vec=2*pi*(0:TrainParam.Mtrain*2-1)/(TrainParam.Mtrain*2);
    f=Amp*[(a-b)*cos(psi_vec)+h*cos((a-b)/b*psi_vec); (a-b)*sin(psi_vec)-h*sin((a-b)/b*psi_vec)];
    plot(f(1,:),f(2,:),'xr')
    hold off
    pause(.5)
end
%% Plot the low dim manifold using PCA
% Calculate correlations
figure
z=squeeze(phi(Ysave(:,:,size(Ysave,3))));
dz=z-repmat(mean(z,2),1,size(z,2));
Cor=1/size(z,1).*(z* z');
[V D]=eig(Cor);
subplot 211
plot(cumsum(flipud(diag(D)))./sum(flipud(diag(D))))
xlim([0 10])
title('Eigenvalues of the correaltion matrix')
subplot 212
plot3(V(:,end)'*dz,V(:,end-1)'*dz,V(:,end-2)'*dz,'o')
title('1st and 2nd PCs')
%% Plot the final state for several angles

disc=0:N-1; theta=2*pi/N*disc;
% TestParam.Mtest=30;
psi_vec_test=2*pi*(0:TestParam.Mtest-1)/TestParam.Mtest;psi_ic=psi_vec_test;
Ang=8;
rateOL=1;
IC=rateOL+Amp*(AmpC*cos(repmat(theta.',1,length(psi_ic))).*cos(repmat(psi_ic,length(theta),1))...
    +AmpS*sin(repmat(theta.',1,length(psi_ic))).*sin(repmat(psi_ic,length(theta),1)));

% IC=rateOL+AmpC*cos(repmat(theta.',1,length(psi_ic))-(repmat(psi_ic,length(theta),1)))  );
figure(55)
for ttt=size(Ysave,3)
    AngInd=0;
    for Ang=1:3:TestParam.Mtest
        AngInd=AngInd+1;
        subplot (5,2,AngInd)
        plot(squeeze(phi(Ysave(:,Ang,ttt)) ))
        hold on
        plot(IC(:,Ang),'r')
        plot(0.*ones(length(IC),1),'--k')
        title(sprintf('time=%d',ttt))
        hold off
        pause(.5)
        
    end
end

