clear;close all
saveFlag=0;
result=struct;
tmax=1000;
dt=1;
rng(1)

g_vec=[0.5 1 1.5];
result.with_learning={};
for TT=1:10
    TT
    tic;
    for GG=1:length(g_vec)
        
        net.g=g_vec(GG);
        net.g
        tic;
        net.N=1000;
        net.phi=@(x)erf(x./sqrt(2));
        net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
%         hp.sim_resolution=120;
        hp.sim_resolution=120;


        [hp, net, sim] = prep_network_param(hp, net,struct);
        % net = update_net_g(net,net.g); % this is a source for bugs and you generate two random matrices.
        %
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        sim.r = net.phi(x);
        size(sim.r )
        
        M_vec=[3,6,8,10,12,15, 20]; % sim_resolution must be a multiplier of 2M
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
                'dt',dt));
            
            
        end
    end
            toc;
end
%%
for GG=1:length(g_vec)
    for MM=1:length(M_vec)
        for TT=1:10
            vel2(GG,MM,TT)=mean(mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:2)),2*pi),[],2),2).^2);
            
            vel2Max(GG,MM,TT)=max(mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:2)),2*pi),[],2),2).^2);
            
            %     vel{MM}=mean(diff((angle(result.with_learning{MM}.z_cl(:,1:2))),[],2),2);% initial psi dot
            %     vel{MM}=mean(diff((angle(result.with_learning{MM}.z_cl(:,1:end))),[],2),2); % integrate psi dot
        end
    end
end
% Plot Drift Vs g,M
%%
figure(5)
% subplot(2,1,1)
for gg=1:3
%     subplot(2,1,1)
    
    % errorbar(g_vec,mean(sqrt(vel2(:,mm,:)),3),std(sqrt(vel2(:,mm,:)),[],3),'o-')
%     semilogy(2*M_vec,squeeze(sqrt(vel2(gg,:,:))),'o-')
        errorbar(2*M_vec,mean(sqrt(vel2(gg,:,:)),3),std(sqrt(vel2(gg,:,:)),[],3),'o-')
set(gca,'yscale','log')
    ylabel('average drift velocity')
    hold all
    box off
%     subplot(2,1,2)
%     % errorbar(g_vec,mean(sqrt(vel2Max(:,mm,:)),3),std(sqrt(vel2Max(:,mm,:)),[],3),'o-')
%     loglog(M_vec,squeeze(sqrt(vel2Max(gg,:,:))),'o-')
    
    xlabel('M')
%     ylabel('Max drift velocity')
    hold all
    box off
    legend('g=0.5','g=1','g=1.5')
end
xlim([0 40])
%% Compare differences between reconstruction error and drift
TT=1;
for MM=1:length(M_vec)
    for GG=1:length(g_vec)
        
        %         v=mean(diff(mod(angle(result.with_learning{GG,MM}.z_cl(2:end,1:2)),2*pi),[],2),2);
        v=mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:2)),2*pi),[],2),2);
        
        v100=mean(diff(mod(angle(result.with_learning{GG,MM,TT}.z_cl(2:end,1:100)),2*pi),[],2),2);
        
        Er=phase([1,1i]*result.with_learning{GG,MM,TT}.z_ol)-sim.psi;
        Er(1)=[];
        overlap(GG,MM)=(Er*v)./(sqrt(Er*Er')*sqrt(v'*v));
        overlap100(GG,MM)=(Er*v100)./(sqrt(Er*Er')*sqrt(v100'*v100));
        Scaling(GG,MM)=norm(v)./norm(Er);
%         Scaling(GG,MM)=mean(v./Er');

%                 plot(Er,'b')
%                 hold on
%                 plot(v,'r')
%                 title(sprintf('M=%d, g=%.1f, overlap=%.3f',2*M_vec(MM),g_vec(GG),overlap(GG,MM)))
%                 pause
%                 hold off
    end
end

figure(10)
subplot 211
semilogy(g_vec,overlap,'o-')
box off
ylim([0.99 1])
legend('M=6','M=12','M=40')
title('Dot product of Reconstruction error and initial velocity')
subplot 212
semilogy(g_vec,overlap100,'o-')
box off
% ylim([0.9 1])
legend('M=6','M=12','M=40')
title('Dot product of Reconstruction error and averaged 100 tau velocity')
xlabel('g')

figure(101)
plot(g_vec,Scaling,'o-')
hold on
box off
% ylim([0.99 1])
legend('M=6','M=12','M=40')
title('Scaling of Reconstruction error and initial velocity')
load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/CompareMF/sweep_leading_pole_vs_g.mat')
plot(J_vec,uu,'--k')

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

