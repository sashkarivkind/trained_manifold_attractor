clear
cnt=0;
load pole_data_vs_g_v1.mat
psi_report_vec=[45]%,30,45,60];
figure;
for psipsi=1:length(psi_report_vec)
    % figure;
            acac=[];

    psi_report=psi_report_vec(psipsi)*pi/180;
tau_eff=[];
    seeds=1:10;
    
    for seed_index=1:length(seeds)
        iter=load(['resultsN4000seed',num2str(seeds(seed_index))]);
        fol=[1,1i]*iter.sim.f_ol;
        fol_new= fol+iter.epsilonmag*(cos(iter.psi)+1i*sin(iter.psi));
        tt0=ceil(iter.tswitch/iter.dt);
        dt=iter.dt;
        zdelta_rec=iter.zdelta_rec;
        gg_to_plot=1:1:(length(iter.g_vec));

        % gg_to_plot=[0 0.5 1 1.5];
        vv=[];
        for gg=gg_to_plot
           vv_vec=diff(angle(zdelta_rec{gg}))/dt;
            [acac(seed_index,gg),ii]=min(abs(angle(zdelta_rec{gg})-psi_report)); %picking index closest to psi_report
            vv(end+1)=vv_vec(min(ii,length(vv_vec))); %adding velocity at that index
            
            cnt=cnt+1
        end
            plot(iter.g_vec,vv/vv(1),'o') % assuming that all the iterations have the same g_vector
    hold on;
    tau_eff(seed_index,:)=1./(vv/vv(1));
    end
   %%

end
plot([0,J_vec],[1,uu],'-');
%%
figure;
plot([0,J_vec],1./[1,uu],'-');
hold on;
tau_eff_mean=mean(tau_eff,1);
tau_eff_stde=std(tau_eff,[],1)/sqrt(size(tau_eff,1)-1);

errorbar(iter.g_vec,tau_eff_mean,tau_eff_stde,'o')
max(abs(acac(:)))