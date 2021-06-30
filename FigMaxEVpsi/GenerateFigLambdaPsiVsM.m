
% Plot lambda max Vs M for different g, eig Vs MF theory
load('maxEv_10_Mpt_sim_interations_N_8000_M_4_6_8_10_12_FAST_diffg.mat')
figure(2)
M_vec=[2 3 4 5 6 7];
    set(gca,'ColorOrderIndex',1);
    zMaxEV=[];
for gg=1:length(g_vec)
    for MM=1:size(results.maxev,1)
        zMaxEV{MM}=((squeeze(results.maxev(MM,gg,1:M_vec(MM),:))));
        mMaxEV(MM)=mean(zMaxEV{MM}(:));
        sdMaxEV(MM)=std(zMaxEV{MM}(:));
        semMaxEV(MM)=std(zMaxEV{MM}(:))./sqrt(length(zMaxEV{MM}(:)));
    end
    semilogy(2*M_vec,abs(mMaxEV),'o');
        box off; xlabel('M');ylabel('\lambda_{max}')
        box off; xlabel('M');ylabel('|\lambda_{max}|')
    hold on
end

legend('g=0.2','g=0.5','g=1','g=1.6','g=1.8')
ylim([0.00001 0.1])
xlim([3.5 14.5])

samplingMFT_overM_with_load
