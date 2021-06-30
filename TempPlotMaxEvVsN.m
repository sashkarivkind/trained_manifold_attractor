% load('maxEv_10_Mpt_sim_interations_N_16k_g_1_8_M_6_8_10_12.mat')
% g_vec=[1.8];
% M_vec=[3 4 5 6];
% figure
% MM=0;
% for M=M_vec
%     MM=MM+1;
% zMaxEV=squeeze(results.maxev(MM,:,1:M,:));
% mMaxEv=mean(zMaxEV(:));
% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% % plot(g_vec,mMaxEv,'-x')
% hold on
% end

%% load 1k
% load('maxEv_10_Mpt_sim_interations_N_1k_g_1_8_M_6_8_10_12.mat')
load('maxEv_10_Mpt_sim_interations_N_1000_M_4_6_8_10_12_FAST.mat')
g_vec=[1.8];
M_vec=[2 3 4 5 6];
% figure
MM=0;
N=1
for M=M_vec
    MM=MM+1;
zMaxEV=squeeze(results.maxev(MM,:,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));
% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% plot(g_vec,mMaxEv,'-x')
end
%% load 2k
load('maxEv_10_Mpt_sim_interations_N_2000_M_4_6_8_10_12_FAST.mat')
g_vec=[1.8];
M_vec=[2 3 4 5 6];
% figure
MM=0;
N=2
for M=M_vec
    MM=MM+1;
zMaxEV=squeeze(results.maxev(MM,:,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));
% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% plot(g_vec,mMaxEv,'-x')
end
%% load 4k
load('maxEv_10_Mpt_sim_interations_N_4000_M_6_12.mat')
g_vec=[1.8];
M_vec=[2 3 4 5 6];
% figure
N=3
MM=2;M=3;
zMaxEV=squeeze(results.maxev(1,end-1,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));
MM=5;M=6;
zMaxEV=squeeze(results.maxev(2,end-1,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));

load('maxEv_10_Mpt_sim_interations_N_4000_M_8_10.mat')
MM=3;M=4;
zMaxEV=squeeze(results.maxev(1,end-1,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));
MM=4;M=5;
zMaxEV=squeeze(results.maxev(2,end-1,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));

% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% plot(g_vec,mMaxEv,'-x')
% end
%% load 8k
load('maxEv_10_Mpt_sim_interations_N_8000_M_4_6_8_10_12_FAST.mat')
g_vec=[1.8];
M_vec=[2 3 4 5 6];
% figure
MM=0;
N=4
for M=M_vec
    MM=MM+1;
zMaxEV=squeeze(results.maxev(MM,:,1:M,:));
mMaxEv(N,MM)=mean(zMaxEV(:));
% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% plot(g_vec,mMaxEv,'-x')
end
%% Load 16k
load('maxEv_10_Mpt_sim_interations_N_16k_g_1_8_M_6_8_10_12.mat')
g_vec=[1.8];
M_vec=[ 3 4 5 6];
% figure
MM=0;
N=5
for M=M_vec
    MM=MM+1;
zMaxEV=squeeze(results.maxev(MM,:,1:M,:));
mMaxEv(N,MM+1)=mean(zMaxEV(:));
% sdMaxEv=std(zMaxEV(:));
% errorbar(g_vec,mMaxEv,sdMaxEv./sqrt(10*M),'-x')
% plot(g_vec,mMaxEv,'-x')
end
%%
mMaxEv(mMaxEv==0)=nan;
figure; plot([1000 2000 4000 8000 16000], mMaxEv,'-o')
legend('M=4','M=6','M=8','M=10','M=12')
box off
title('g=1.8; N=1k,2k,8k- eigs')
 xlabel('N'); ylabel('\lambda_{max}')