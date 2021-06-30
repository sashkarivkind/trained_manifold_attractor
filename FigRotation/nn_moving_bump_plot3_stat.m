figure;
tt0=ceil(tswitch/dt);
cnt=1;
gg_to_plot=1:length(g_vec);%[1,3,6,8];
    for gg=gg_to_plot
        %     vv(end+1)=max(diff(angle(zdelta_rec{gg})))/dt;
        vv_vec=diff(angle(zdelta_rec{gg}))/dt;
%         [~,ii]=min(abs(angle(zdelta_rec{gg})-psi_report));
%         vv(end+1)=vv_vec(min(ii,length(vv_vec)));
        v_vec=diff(angle(zdelta_rec{gg}(tt0:end)))/dt;
        t_vec=dt*(0:(length(zdelta_rec{gg})-1));
            h(cnt)=plot(t_vec, angle(zdelta_rec{gg})/pi*180,'linewidth',2);
        % %     plot(angle(zdelta_rec{gg}));
            hold on;
%             set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
%             plot((0:(length(fol)-1))*360/160,...
%                 uu1g(gg)*(angle(fol_new)-angle(fol)),':','linewidth',2);
    cnt=cnt+1;
    end
    
%     plot(g_vec,vv/vv(1),'*')
%     hold on;
% end
% plot(J_vec,uu,'-');
% prep_legend = ...
%     arrayfun(@(x) [num2str(x),'^{\circ}'],...
%     90-psi_report_vec','UniformOutput',0)
% legend([prep_legend; 'theory']);
% xlim([0,90]);
box off;
legend(h,num2str(g_vec(gg_to_plot)'));