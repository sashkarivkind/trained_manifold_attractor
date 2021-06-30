% clear;
tmax=1000;

N=4000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);

g_vec=[0.01, 0.1, 0.3, 0.5];
A_vec=[0.5,1,1.2,1.5,2.];
%% cover up for forgetting to save J2
uu = load('A_J2_table');
J2_vec=[];
for aa = 1:length(A_vec)
 ii=find(uu.A_vec==A_vec(aa));
 if length(ii)~=1
     error
 end
 J2_vec(ii)=uu.J2_vec(ii);
end
%%

seed_vec=1:5;
rmse_rec=[];
A_rec=[];
g_rec=[];
seed_rec=[];
rmsv_rec=[];
cnt=0;
for aa=1:length(A_vec)
    A=A_vec(aa);
    for gg=1:length(g_vec)
        g=g_vec(gg);
        for sese=1:length(seed_vec)
            try
            loaded=load(['results/ZZZresultsN',num2str(N),'iter',num2str(cnt)],'result','N','g','A','seed');
            cnt=cnt+1;
            rmse_rec(end+1)=loaded.result.rmse;
            A_rec(end+1)=loaded.A;
            g_rec(end+1)=loaded.g;
            seed_rec(end+1)=loaded.seed;
            rmsv_rec(end+1)=sqrt(mean((loaded.result.v.^2)));
            catch
                break
            end
        end
    end
end
%%
% figure; 
% dmat=cellfun(@(xx) xx.d, mft_result_table);
% Amat=cellfun(@(xx) xx.A, mft_result_table);
% rmsmat=dmat./Amat/sqrt(N);
% for aa=3%1:length(A_vec)
% cc1=abs(A_rec-A_vec(aa))<1e-8;
% plot(g_rec(cc1),rmse_rec(cc1),'x');
% hold on;
% set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
% plot(g_vec,rmsmat(aa,:),'-');
% end

%%
figure;
subplot(2,2,1)
dmat=cellfun(@(xx) xx.d, mft_result_table);
Amat=cellfun(@(xx) xx.A, mft_result_table);
rmsmat=dmat./Amat/sqrt(N);
for gg=1:length(g_vec)
    subplot(2,2,gg)
cc1=abs(g_rec-g_vec(gg))<1e-8;
plot(A_rec(cc1),rmse_rec(cc1),'x');
hold on;
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
plot(A_vec,rmsmat(:,gg),'-');
xlabel('A_0');
ylabel('\Delta_\Psi');
title(['g=',num2str(g_vec(gg))]);
end
% %%
% figure; 
% for aa=1:length(seed_vec)
% cc1=abs(seed_rec-seed_vec(aa))<1e-8;
% plot(g_rec(cc1),rmse_rec(cc1),'x');
% hold on;
% end
% %%
%%
% figure;
% plot(rmse_rec,rmsv_rec,'x');
% hold on;
% plot([0,1],[0,1],'k--');
% xlim([0,max([rmse_rec,rmsv_rec])]);
% ylim([0,max([rmse_rec,rmsv_rec])]);
% %%
% figure;
% loglog(rmse_rec,rmsv_rec,'x');
% hold on;
% loglog([1e-10,1],[1e-10,1],'k--');
% xlim([min([rmse_rec,rmsv_rec]),max([rmse_rec,rmsv_rec])]);
% ylim([min([rmse_rec,rmsv_rec]),max([rmse_rec,rmsv_rec])]);