%%
% load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/fig_prep_scripts/Gp_10_sim_interations_N_4000.mat')
% load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/fig_prep_scripts/Gp_10_sim_interations_N_4000_M_12')
load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/fig_prep_scripts/Gp_10_sim_interations_N_16000_M_6_Full.mat');

% load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/fig_prep_scripts/Gp_10_sim_interations_N_16000_M_12_Full.mat');

tmax=1000;
dt=1;
g_vec=[0.2:0.2:2];
% g_vec=[1.2:0.2:2];
% M_vec=[3 5];
 M_vec=[6];
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;

% figure;
G_vec_over_m_g={};
for MM=1
    for GG=1:length(g_vec)
        for II=1:10
            for pp=1:length(results.G_polar{MM,GG,II}.flat)
                G_vec_over_m_g{MM}(GG,pp,II)=results.G_polar{MM,GG,II}.flat{pp}{2,2};
            end
        end
    end
end
%
% plot(g_vec,squeeze(G_vec_over_m_g{1}(:,:)),'.r')
hold all
errorbar(g_vec,mean(G_vec_over_m_g{1}(:,:),2),std(G_vec_over_m_g{1}(:,:),[],2)./sqrt(10),'or','MarkerSize',3)
hold on;
% plot(g_vec,squeeze(G_vec_over_m_g{2}(:,:)),'.b')
% errorbar(g_vec,mean(G_vec_over_m_g{2}(:,:),2),std(G_vec_over_m_g{2}(:,:),[],2)./sqrt(10),'ob')

load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/mft_fun/Gp_MF_M6_M10.mat')
% load('/Users/darshanr/Dropbox (HHMI)/LearningForce 2017/BalanceLearnSasha/nn2021/mft_fun/Gp_MF_M12.mat')

plot(g_vec,G1p(1,:),'-r')
% plot(g_vec,G1p(2,:),'-b')

plot(0:0.1:2,ones(length(0:0.1:2),1),'-k')

ylim([0.8 1.05])
xlim([0 2.1])
box off
axis square