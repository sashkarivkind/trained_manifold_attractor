clear;
% close all;
saveFlag=0;
RunCLFlag=1;

saveFolder='/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/FiguresSubPlots/FigLearnExample/';
result=struct;
tmax=1000;
dt=1;
UniformFlag=0;
% a=2;b=1;A=1.5;rot=0;h=0.5;
a=4;b=1;A=1.2;rot=0;h=0.0;
g_vec=[1.5];

net.g=0; net.N=1000;
net.phi=@(x)erf(x./sqrt(2)); net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp=struct;
sim=struct;
hp.A=A;

rng(1)

trained_net={};

[hp, net, sim] = prep_network_param(hp, net,sim);

if UniformFlag
    tmp_f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
    sim.f_ol=curvspace(tmp_f_ol',hp.sim_resolution);
    sim.f_ol=sim.f_ol';
else
    sim.f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
end

% g_vec=[0.01 0.1:0.1:1.4];
x0 = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));



for gg=1:length(g_vec)
    g=g_vec(gg)
    trained_net{gg} = update_net_g(net,g);

    net
    hp
    sim
    x = fast_conv_to_fp(trained_net{gg},sim.f_ol,struct('ol',1));
    sim.r = trained_net{gg}.phi(x);
    pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
%     trained_net{gg}.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
               regfac = hp.M*hp.alpha_reg*eye(length(pts));
                trained_net{gg}.wout = (sim.r(:,pts)/...
                    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
                
                
    result.with_learning{gg}.z_ol = trained_net{gg}.wout'*sim.r;
    
        [result.with_learning{gg}.x_cl,result.with_learning{gg}.z_cl]=fast_conv_to_fp_extended(trained_net{gg},[],...
            struct('xinit',x,...
            'ol',0,...
            'tmax',tmax,...
            'dt',dt));
    % %
    % %
    % %     if RunCLFlag
    %         %simulating closed loop
    x_test_rand=fast_conv_to_fp(trained_net{gg},[],struct('xinit',5*randn(size(x))));
    r_test_rand=trained_net{gg}.phi(x_test_rand);
    sim.z_test_rand= trained_net{gg}.wout'*r_test_rand;
    
    x_test_noi=fast_conv_to_fp(trained_net{gg},[],struct('xinit',x+1e-3*randn(size(x))));
    r_test_noi=trained_net{gg}.phi(x_test_noi);
    sim.z_test_noi= trained_net{gg}.wout'*r_test_noi;
    
    %plotting results
    
    plot(result.with_learning{gg}.z_ol(1,:),result.with_learning{gg}.z_ol(2,:),'.');
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
    plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
    plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
    plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
    
    xlim([-2 2]); ylim([-2 2])
    %     end
end

% figure(33)
id=[1:50:1000];
% id=1:1000;
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
% figMemTuning2=figure(34);
figure
% plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(end/2+id,:)))');
plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(id,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])


fuRand=fft(net.phi(result.with_learning{gg}.x_cl(id,:)),[],2);

% f1_stat_id=abs(fuRand(:,2))./size(fuRand,2)
%     title(['g=1.4 f_1=',num2str(f1_stat_id')])

title(['a= ',num2str(a), ', h= ',num2str(h), ', g= ',num2str(g) ])
% [m mm]=sort(f1_stat_id);
%% Calculate statistics
close all

fuBase=fft(net.phi(result.baseline.x_cl(1:end,:)),[],2);
figure
for gg=1:length(g_vec)
fuRand=fft(net.phi(result.with_learning{gg}.x_cl(1:end,:)),[],2);

f1_stat(gg,:)=abs(fuRand(:,2))./size(fuRand,2);
sdf1_stat(gg)=std(f1_stat(gg,:));
subplot (4,4,gg)
hist(f1_stat(gg,:),50)
xlim([0 1])
title(num2str(g_vec(gg)))
xlabel('f_1'); ylabel('counts')
box off

end
figure(100)
plot(g_vec,sdf1_stat,'ob');hold on
box off
xlabel('heteroginity level (g)'); ylabel ('Tuning curves diversity (SD (f_1) )');
%% Generate MF estimate of te sd(f^1_i)
N=1000;
tpN=[0:(N-1)]/N*2*pi;tpN=tpN';
Amp=hp.A; res=160;
psi_vec=[0:(res-1)]/res*2*pi;
f1_stat_MF=[];
sdf1_stat_MF=[];
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
    result.with_learning{gg}.MF_TC=[];
    for ii=1:length(psi_vec)
        harmonics=psi_vec(ii)'*[1:2:length(c)];
        x1_approx=sum((ai.*cos(harmonics))...
            +(bi.*sin(harmonics)),2);
        result.with_learning{gg}.MF_TC(:,ii)=net.phi(Amp*cos(tpN-psi_vec(ii))+x1_approx);
    end

    % Calculate f_1
    fuRandMF=fft(result.with_learning{gg}.MF_TC,[],2);
    f1_stat_MF(gg,:)=abs(fuRandMF(:,6))./size(fuRandMF,2);
    sdf1_stat_MF(gg)=std(f1_stat_MF(gg,:));
    % Plot the statistics
    figure(33)
    subplot (4,4,gg)
    hist(f1_stat_MF(gg,:),50)
    xlim([0 1])
    title(num2str(g_vec(gg)))
    xlabel('f_1'); ylabel('counts')
    box off

end

figure(100)
plot(g_vec,sdf1_stat_MF,'-k');hold on
% plot(g_vec,sdf1_stat_MF_test,'-r');hold on
%Plot Tuning curves example with high and low SI

gg=1 %0.2;
% id=[20 21 22];
% id=[520 182 521 609];
%%
figure(33)
id=[1:150:1000];
% id=1:1000;
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
figMemTuning2=figure(24);
% plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(end/2+id,:)))');
plot(xaxis_deg-pi,(net.phi(result.with_learning{gg}.x_cl(id,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])


fuRand=fft(net.phi(result.with_learning{gg}.x_cl(id,:)),[],2);

% f1_stat_id=abs(fuRand(:,2))./size(fuRand,2)
%     title(['g=1.4 f_1=',num2str(f1_stat_id')])

title(['a= ',num2str(a), ', h= ',num2str(h), ', g= ',num2str(g) ])
% [m mm]=sort(f1_stat_id);
%% Calculate participation ratio

