%need to run nn_moving_bump_plot1 first to make this script work correctly
clear;
% figure;
cnt=1;

load('resultsN4000seed2','g_vec');
% uu1g=[1, uu(3),mean(uu(5:6)),uu(8)];

%%
ol_pole_vec=zeros(size(g_vec));
beta_info=load('../mft_fun/data/beta_graph_long.mat');
for gg=1:length(g_vec)
    angle_data{gg}=[];
    v_data{gg}=[];
    ii=find(abs(beta_info.J_vec-g_vec(gg))<0.02);
    if length(ii)~=1
        error('pole info not found!');
    end
    ol_pole_vec(gg)=beta_info.debuu{ii}.pp_p1.v1(1);
end
%%

gg_to_plot=1:1:(length(g_vec));
seeds=1:10;

for seed_index=1:length(seeds)
    iter=load(['resultsN4000seed',num2str(seeds(seed_index))]);
    fol=[1,1i]*iter.sim.f_ol;
    fol_new= fol+iter.epsilonmag*(cos(iter.psi)+1i*sin(iter.psi));
    tt0=ceil(iter.tswitch/iter.dt);
    dt=iter.dt;
    zdelta_rec=iter.zdelta_rec;
    % gg_to_plot=[0 0.5 1 1.5];
    for gg=gg_to_plot
        v_vec=diff(angle(zdelta_rec{gg}(tt0:end)))/dt;
        angle_data{gg}(seed_index,:)=angle(zdelta_rec{gg}(tt0+1:end))/pi*180;
        v_data{gg}(seed_index,:)=v_vec;
        %     h(cnt)=plot(angle(zdelta_rec{gg}(tt0+1:end))/pi*180,v_vec,'linewidth',2);
        %     hold on;
        %     set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
        %     plot((0:(length(fol)-1))*360/160,...
        %         ol_pole_vec(gg)*(angle(fol_new)-angle(fol)),':','linewidth',2);
        cnt=cnt+1
    end
end
%%
figure;
cnt=0;
h=[];
    for gg=1:5:21;%gg_to_plot
cnt=cnt+1
%            h(cnt)=errorbar(angle_data{gg}(1,:),...
%                mean(v_data{gg}),std(v_data{gg}),'linewidth',2);
ccc=get(gca,'ColorOrderIndex');
%            h(cnt)=errorbar(mean(angle_data{gg}),...
%                mean(v_data{gg}),std(v_data{gg}),'linewidth',2);


      h(cnt)=plot(mean(angle_data{gg}),...
               mean(v_data{gg}),'linewidth',2);

%            plot(angle_data{gg}',v_data{gg}','linewidth',2);
            hold on;
            set(gca,'ColorOrderIndex',ccc);
            plot((0:(length(fol)-1))*360/160,...
                ol_pole_vec(gg)*(angle(fol_new)-angle(fol)),':','linewidth',2);
    end
    %%
        for gg=gg_to_plot
        %     vv(end+1)=max(diff(angle(zdelta_rec{gg})))/dt;
%         [~,ii]=min(abs(angle(zdelta_rec{gg})-psi_report));
%         vv(end+1)=vv_vec(min(ii,length(vv_vec)));
        t_vec=dt*(0:(length(zdelta_rec{gg})-1));
            h(cnt)=plot(t_vec, angle(zdelta_rec{gg})/pi*180,'linewidth',2);
        % %     plot(angle(zdelta_rec{gg}));
            hold on;
%             set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
%             plot((0:(length(fol)-1))*360/160,...
%                 uu1g(gg)*(angle(fol_new)-angle(fol)),':','linewidth',2);
    cnt=cnt+1;
    end
    
% box off;
% legend(h,num2str(g_vec(gg_to_plot)'));