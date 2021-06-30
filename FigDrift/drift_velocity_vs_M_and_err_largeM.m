clear;close all
saveFlag=0;
result=struct;
tmax=3;
dt=.5;
rng(1)

% g_vec=[0.5 1 1.5];
g_vec=[0.2 0.5 1.0 1.6 1.8];
%         M_vec=[3,6,8,10,12,15, 20]; % sim_resolution must be a multiplier of 2M
% M_vec=4:2:14; % sim_resolution must be a multiplier of 2M
M_vec=[2:6 10 20 30 40 60];%2:7; % sim_resolution must be a multiplier of 2M

M_legend=eq_legend(2*M_vec,'M=');
g_legend=eq_legend(g_vec,'g=');
%%
ol_pole_vec=zeros(size(g_vec));
beta_info=load('../mft_fun/data/beta_graph_long.mat');
for gg=1:length(g_vec)
    ii=find(abs(beta_info.J_vec-g_vec(gg))<1e-10);
    if length(ii)~=1
        error('pole info not found!');
    end
    ol_pole_vec(gg)=beta_info.debuu{ii}.pp_p1.v1(1);
end
%%
% s
result.with_learning={};

for TT=1:10
    TT
    tic;
    for GG=1:length(g_vec)
        
        net.g=g_vec(GG);

        tic;
        net.N=4000;
        net.phi=@(x)erf(x./sqrt(2));
        net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
%         hp.sim_resolution=120;
        hp.sim_resolution=240;


        [hp, net, sim] = prep_network_param(hp, net,struct);
        % net = update_net_g(net,net.g); % this is a source for bugs and you generate two random matrices.
        %
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        sim.r = net.phi(x);
        size(sim.r )
        

        for MM=1:length(M_vec)
            
            hp.M=M_vec(MM);
            
            pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
            regfac = hp.M*hp.alpha_reg*eye(length(pts));
            net.wout = (sim.r(:,pts)/...
                (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
            
            result.with_learning{GG,MM,TT}.z_ol = net.wout'*sim.r;
            [x_cl,result.with_learning{GG,MM,TT}.z_cl]=fast_conv_to_fp_extended(net,[],...
                struct('xinit',x,...
                'ol',0,...
                'tmax',tmax,...
                'dt',dt,'save_neurons',1));
            
            
        end
    end
            toc;
end
%%
for GG=1:length(g_vec)
    for MM=1:length(M_vec)
        for TT=1:10
            vel2(GG,MM,TT)=1/dt^2*mean(mean(...
                diff( mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,:)),2*pi),[],2).^2));
            
            vel2Max(GG,MM,TT)=1/dt^2*max(mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,:)),2*pi),[],2),2).^2);
            
%             vel2(GG,MM,TT)=1/dt^2*mean(...
%                 diff( mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,[1 11])),2*pi),[],2).^2);
%             vel2Max(GG,MM,TT)=1/dt^2*max(mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,[1 11])),2*pi),[],2),2).^2);
            
            %     vel{MM}=mean(diff((angle(result.with_learning{MM}.z_cl(:,1:2))),[],2),2);% initial psi dot
            %     vel{untitled.figMM}=mean(diff((angle(result.with_learning{MM}.z_cl(:,1:end))),[],2),2); % integrate psi dot
        end
    end
end
%% Compare differences between reconstruction error and drift
for TT=1:10
for MM=1:length(M_vec)
    for GG=1:length(g_vec)
        
        %         v=mean(diff(mod(angle(result.with_learning{GG,MM}.z_cl(2:end,1:2)),2*pi),[],2),2);
        v=1/dt*mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:2)),2*pi),[],2),2);
        
%         v100=mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:100)),2*pi),[],2),2);
        
        Er=phase([1,1i]*result.with_learning{GG,MM,TT}.z_ol)-sim.psi;
        Errec(GG,MM,TT)=sqrt(mean(Er.^2));
        Er(1)=[];
        overlap(GG,MM)=(Er*v)./(sqrt(Er*Er')*sqrt(v'*v));
%         overlap100(GG,MM)=(Er*v100)./(sqrt(Er*Er')*sqrt(v100'*v100));
        Scaling(GG,MM)=norm(v)./norm(Er)
 
%         Scaling(GG,MM)=mean(v./Er');

%                 plot(Er,'b')
%                 hold on
%                 plot(v,'r')
%                 title(sprintf('M=%d, g=%.1f, overlap=%.3f',2*M_vec(MM),g_vec(GG),overlap(GG,MM)))
%                 pause
%                 hold off
    end
end
end
%%
% load('../../../BalanceLearnSasha/CompareMF/sweep_leading_pole_vs_g.mat','J_vec','uu')
% % plot(J_vec,uu,'--k')
% clear slow_down_fac
% for gg=1:length(g_vec)
% slow_down_fac(gg)=uu(abs(J_vec-g_vec(gg))<1e-10);
% end
%% Plot Drift Vs g,M
%
figure
hh=[];
for gg=1:length(g_vec)
        hh(end+1)=errorbar(2*M_vec,mean(sqrt(vel2(gg,:,:)),3),std(sqrt(vel2(gg,:,:)),[],3),'o-');
        hold on;     
        set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
        plot(2*M_vec,ol_pole_vec(gg)*Errec(gg,:),'*--')
set(gca,'yscale','log')
    ylabel('average drift velocity')
    hold all
    box off  
    xlabel('M')
    hold all
    box off
    legend(hh,g_legend)
    
%     figure(6)
%      errorbar(2*M_vec,mean(sqrt(vel2Max(gg,:,:)),3),std(sqrt(vel2Max(gg,:,:)),[],3),'o-')
% set(gca,'yscale','log')
%     hold all
%     box off
% 
%     xlabel('M')
%     ylabel('Max drift velocity')
%     hold all
%     box off
%     legend(g_legend)
end
% xlim([0 40])
%%
figure(10)
subplot 211
semilogy(g_vec,overlap,'o-')
box off
ylim([0.99 1])
legend(M_legend)
title('Dot product of Reconstruction error and initial velocity')
subplot 212
% semilogy(g_vec,overlap100,'o-')
box off
% ylim([0.9 1])
legend(M_legend)
title('Dot product of Reconstruction error and averaged 100 tau velocity')
xlabel('g')

figure(101)
plot(g_vec,Scaling,'o-')
hold on
box off
% ylim([0.99 1])
legend(M_legend)
title('Scaling of Reconstruction error and initial velocity')
% load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/CompareMF/sweep_leading_pole_vs_g.mat')


%%
Er=phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on;
plot(sim.psi-pi,vel{1},'--r');
M=M_vec(1)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')


box off;
%ylim([-0.03 0.03]);xlim([-pi pi])



if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM6g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM6g0_5rnd1N2k.pdf')
    end
end

Er=phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on;
plot(sim.psi-pi,vel{2},'--r');
M=M_vec(2)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')
box off;%ylim([-0.03 0.03]);
xlim([-pi pi])


if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM12g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM12g0_5rnd1N2k.pdf')
    end
end

Er=phase([1,1i]*result.with_learning{3}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on;
plot(sim.psi-pi,vel{3},'--r');
M=M_vec(3)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')
box off;%ylim([-0.03 0.03]);
xlim([-pi pi])


if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM40g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM40g0_5rnd1N2k.pdf')
    end
end

figure;plot(angle(result.with_learning{1}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])
if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM6g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM6g0_5rnd1N2k.pdf')
    end
end
figure;plot(angle(result.with_learning{2}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])
if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM12g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM12g0_5rnd1N2k.pdf')
    end
end
figure;plot(angle(result.with_learning{3}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])

if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM40g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM40g0_5rnd1N2k.pdf')
    end
end

