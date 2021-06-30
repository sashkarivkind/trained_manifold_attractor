%need to run nn_moving_bump_plot1 first to make this script work correctly
figure;
tt0=ceil(tswitch/dt);
cnt=1;

fol=[1,1i]*sim.f_ol;
fol_new= fol+epsilonmag*(cos(psi)+1i*sin(psi));
% uu1g=[1, uu(3),mean(uu(5:6)),uu(8)];

%%
ol_pole_vec=zeros(size(g_vec));
beta_info=load('../mft_fun/data/beta_graph_long.mat');
for gg=1:length(g_vec)
    ii=find(abs(beta_info.J_vec-g_vec(gg))<0.02);
    if length(ii)~=1
        error('pole info not found!');
    end
    ol_pole_vec(gg)=beta_info.debuu{ii}.pp_p1.v1(1);
end
%%


% gg_to_plot=1:5:(length(g_vec)-1);

gg_to_plot=[2 6 11 13 15 18];
% gg_to_plot=1:1:(length(g_vec));

% gg_to_plot=[0 0.5 1 1.5];
    for gg=gg_to_plot
        %     vv(end+1)=max(diff(angle(zdelta_rec{gg})))/dt;
        vv_vec=diff(angle(zdelta_rec{gg}))/dt;
%         [~,ii]=min(abs(angle(zdelta_rec{gg})-psi_report));
%         vv(end+1)=vv_vec(min(ii,length(vv_vec)));
        v_vec=diff(angle(zdelta_rec{gg}(tt0:end)))/dt;
            h(cnt)=plot(angle(zdelta_rec{gg}(tt0+1:end))/pi*180,v_vec,'linewidth',2);
%             h(cnt)=plot(angle(zdelta_rec{gg}(tt0+1:end)),v_vec,'linewidth',2);

        % %     plot(angle(zdelta_rec{gg}));
            hold on;
%             set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
            plot((0:(length(fol)-1))*360/160,...
                ol_pole_vec(gg)*(angle(fol_new)-angle(fol)),':','linewidth',2);
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