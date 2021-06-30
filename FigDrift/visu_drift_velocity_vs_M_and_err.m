clear;
close all
% load drift_velocity_vs_M_and_err_RESUME_GAMES
% load drift_velocity_vs_M_and_err_RESUME_GAMES_N4000
load drift_velocity_vs_M_and_err_RESUME_GAMES_N4000_LargeM_10It_B

%% Plot Drift Vs g,M
%
figure
hh=[];
for gg=1:length(g_vec)
%         hh(end+1)=errorbar(2*M_vec,mean(sqrt(vel2(gg,:,:)),3),std(sqrt(vel2(gg,:,:)),[],3),'o-');
        hh(end+1)=plot(2*M_vec,mean(sqrt(vel2(gg,:,1:10)),3),'o');

        hold on;     
        set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
% % % % %                 plot(2*M_vec,ol_pole_vec(gg)*(Errec(gg,:)),'*--')%%%%%
        plot(2*M_vec,ol_pole_vec(gg)*mean(Errec(gg,:,1:10),3),'--')%%%%%
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
