% clear;
% Do the calculations only for the learned angles!
% result=struct;
tmax=1000;
dt=1;
% % g_vec=[0.2:0.2:2];
% g_vec=[0.2 0.5 1 1.6 1.8];
% g_vec=[0.2:0.2:2];
g_vec=[0.2 0.5 1 1.6 1.8];
M_vec=[2 3 4 5 6 7];
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;
results.G_polar={};
%
net.N=8000;

for MM=1:length(M_vec)
    hp.M=M_vec(MM);
    MM
    for GG=1:length(g_vec)
        g=g_vec(GG);
        net.g=g;
        GG
        tic;
        %todo: save or force seed for RNG!
        
        for II=1:10
            result.with_learning={}; Gij={};
             hp.sim_resolution=hp.M;
%             [hp, net, sim] = prep_network_param(hp, net,struct);

            [hp, net, sim] = prep_network_param_learned_angles_only(hp, net,struct);
%             
            
            x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
            sim.r = net.phi(x);
%             pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
            regfac = hp.M*hp.alpha_reg*eye(hp.M);
            net.wout = (sim.r/...
                (sim.r'*sim.r+regfac))*sim.f_ol'; %obtain least mean square solution
            
            cnt=0;
            for pt = 1:hp.M
                cnt=cnt+1;
                S=(net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt)))-eye(net.N);
               results.ev{MM,GG,pt,II}=eigs(S,1,'smallestabs');
               results.maxev(MM,GG,pt,II)=max(real(results.ev{MM,GG,pt,II}));
            end
        end
        toc;
    end
end
%
save('maxEv_10_Mpt_sim_interations_N_8000_M_4_6_8_10_12_FAST_diffg.mat','hp','sim','results','g_vec','M_vec')
%%
% g_vec=[0.2:0.2:2];
% load('maxEv_5_Mpt_sim_interations_N_2000_M_4_6_8_10_12_FAST_diffgAll.mat')
% load('maxEv_10_Mpt_sim_interations_N_4000_M_4_6_8_10_12_FAST_diffg.mat')

load('maxEv_10_Mpt_sim_interations_N_8000_M_4_6_8_10_12_FAST_diffg.mat')
% load('maxEv_1_Mpt_sim_interations_N_16000_M_4_6_8_10_12_FAST_diffg.mat')
figure(2)
M_vec=[2 3 4 5 6 7];
    set(gca,'ColorOrderIndex',1);
    zMaxEV=[];
for gg=1:length(g_vec)
    for MM=1:size(results.maxev,1)
        zMaxEV{MM}=((squeeze(results.maxev(MM,gg,1:M_vec(MM),:))));

%         zMaxEV{MM}=((squeeze(results.maxev(MM,gg,:,:))));
%         zMaxEV{MM}=((squeeze(results.maxev(MM,gg,:,:))));

        mMaxEV(MM)=mean(zMaxEV{MM}(:));
        sdMaxEV(MM)=std(zMaxEV{MM}(:));
        semMaxEV(MM)=std(zMaxEV{MM}(:))./sqrt(length(zMaxEV{MM}(:)));
    end
%     subplot (2,1,1)
%     title('N=4k, eigs')
    
%     errorbar(2*M_vec,mMaxEV,semMaxEV,'o');hold on

    semilogy(2*M_vec,abs(mMaxEV),'o');
    
%     set(gca,'Yscale','log')
    box off; xlabel('M');ylabel('\lambda_{max}')
%     subplot (2,1,2)
%     semilogy(2*M_vec,abs(mMaxEV),'o')
        box off; xlabel('M');ylabel('|\lambda_{max}|')

    hold on
end

legend('g=0.2','g=0.5','g=1','g=1.6','g=1.8')
ylim([0.00001 0.1])
xlim([3.5 14.5])
% set(gca,'yscale','log')
%% N=1000
% g_vec=[0.2:0.2:2];
load('maxEv_5_Mpt_sim_interations_N_2000_M_4_6_8_10_12_FAST_diffgAll.mat')
% load('maxEv_10_Mpt_sim_interations_N_1000_M_4_6_8_10_12_FAST_diffg.mat')
% load('maxEv_10_Mpt_sim_interations_N_4000_M_4_6_8_10_12_FAST_diffg.mat')

figure(2)
M_vec=[2 3 4 5 6 7];
    set(gca,'ColorOrderIndex',1);

for gg=3
    for MM=1:size(results.maxev,1)
        zMaxEV{MM}=((squeeze(results.maxev(MM,gg,1:M_vec(MM),:))));
%         mMaxEV(MM)=mean(zMaxEV{MM}(:));
%         sdMaxEV(MM)=std(zMaxEV{MM}(:));
%         semMaxEV(MM)=std(zMaxEV{MM}(:))./sqrt(length(zMaxEV{MM}(:)));
         semilogy(2*M_vec(MM),abs(zMaxEV{MM}(:)),'og'); hold on
    end
%     subplot (2,1,1)
    title('N=4k, eigs')
    
%     errorbar(2*M_vec,mMaxEV,semMaxEV,'o');hold on

%     semilogy(2*M_vec,abs(mMaxEV),'x');
  
    
%     set(gca,'Yscale','log')
    box off; xlabel('M');ylabel('\lambda_{max}')
%     subplot (2,1,2)
%     semilogy(2*M_vec,abs(mMaxEV),'o')
        box off; xlabel('M');ylabel('|\lambda_{max}|')

    hold on
end

legend('g=0.2','g=0.5','g=1','g=1.6','g=1.8')
% set(gca,'yscale','log')
%%
for MM=1:length(M_vec)
    for II=1:5
        for pt = 1:M_vec(MM)
            results.ev{MM,5,pt,II}=resultsB.ev{MM,1,pt,II};
            results.maxev(MM,5,pt,II)=resultsB.maxev(MM,1,pt,II);
        end
    end
end