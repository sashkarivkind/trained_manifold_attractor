load pole_data_vs_g_v1.mat
psi_report_vec=[1]%,30,45,60];
figure;
for psipsi=1:length(psi_report_vec)
vv=[];
% figure;
psi_report=psi_report_vec(psipsi)*pi/180;
for gg=1:length(g_vec)
%     vv(end+1)=max(diff(angle(zdelta_rec{gg})))/dt;
vv_vec=diff(angle(zdelta_rec{gg}))/dt;
[~,ii]=min(abs(angle(zdelta_rec{gg})-psi_report));
vv(end+1)=vv_vec(min(ii,length(vv_vec)));
    v_vec=diff(angle(zdelta_rec{gg}))/dt;
%     plot(angle(zdelta_rec{gg}(2:end))/pi*180,v_vec);
% %     plot(angle(zdelta_rec{gg}));
%     hold on;
end

 plot(g_vec,vv/vv(1),'o')
hold on; 
end
plot([0,J_vec],[1,uu],'-');
% plot(J_vec,uu2,'-');
prep_legend = ...
    arrayfun(@(x) [num2str(x),'^{\circ}'],...
    psi_report_vec','UniformOutput',0)
legend([prep_legend; 'theory']);
xlim([0,1.6]);
box off;