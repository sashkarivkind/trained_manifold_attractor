% clear;
% Do the calculations only for the learned angles!
% result=struct;
tmax=1000;
dt=1;
g_vec=[0.2:0.2:2];
M_vec=[4 5];
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;
results.G_polar={};
%
net.N=4000;

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
               results.ev{MM,GG,pt,II}=eig(S,'vector');
               results.maxev(MM,GG,pt,II)=max(real(results.ev{MM,GG,pt,II}));
            end
        end
        toc;
    end
end
%
save('maxEv_10_Mpt_sim_interations_N_4000_M_8_10.mat','hp','sim','results')
%%
M_vec=[3 6];
for II=1:10
    for GG=1:10
        for MM=1:2
            M=M_vec(MM);
            for pt = 1:M
                results.maxev(MM,GG,pt,II)=max(real(results.ev{MM,GG,pt,II}));
            end
        end
    end
end
save('maxEv_10_Mpt_sim_interations_N_4000_M_6_12.mat','hp','sim','results')

%%
g_vec=[0.2:0.2:2];
load('maxEv_10_Mpt_sim_interations_N_4000_M_6_12.mat')

zMaxEV=squeeze(results.maxev(1,:,1:3,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(30),'-ob')

hold all
zMaxEV=squeeze(results.maxev(2,:,:,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(60),'-or')

title('10 realizations, average over the M points in each realization')

xlabel('g')
ylabel('\lambda_{max}')

clear
g_vec=[0.2:0.2:2];
load('maxEv_10_Mpt_sim_interations_N_4000_M_8_10.mat')

zMaxEV=squeeze(results.maxev(1,:,1:4,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(40),'-og')

hold all
zMaxEV=squeeze(results.maxev(2,:,:,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(80),'-ok')

legend('M=6','M=12','M=8','M=10')
%% Plot max EV Vs M for different g

g_vec=[0.2:0.2:2];
load('maxEv_10_Mpt_sim_interations_N_4000_M_6_12.mat')
zMaxEV=[];
zMaxEV{2}=squeeze(results.maxev(1,:,1:3,:)); % g, M, it
zMaxEV{2}=squeeze(results.maxev(1,:,1:6,:)); 


mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(30),'-ob')

hold all
zMaxEV=squeeze(results.maxev(2,:,:,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(60),'-or')

title('10 realizations, average over the M points in each realization')

xlabel('g')
ylabel('\lambda_{max}')

clear
g_vec=[0.2:0.2:2];
load('maxEv_10_Mpt_sim_interations_N_4000_M_8_10.mat')

zMaxEV=squeeze(results.maxev(1,:,1:4,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(40),'-og')

hold all
zMaxEV=squeeze(results.maxev(2,:,:,:));
mMaxEv=mean(mean(zMaxEV(:,:,:),2),3);sdMaxEv=std(zMaxEV(:,:),[],2);
errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(80),'-ok')

legend('M=6','M=12','M=8','M=10')



%%

MM=1;
II=8;
for GG=1:length(g_vec)
    plot(results.ev{MM,GG,II},'.')
    hold on 
    plot(results.maxev(MM,GG,II), 0,'ro')
    pause
    hold off
end
