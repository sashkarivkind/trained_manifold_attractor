clear;
tmax=1000;

g_vec=[0.01, 0.1:0.1:2];
A_vec=[0.5,1,1.2,1.5,2.];
seed_vec=1:10;

cnt=0;
%%
% for aa=1:length(A_vec)
for aa=3
    A=A_vec(aa);
    for gg=8:length(g_vec)
        g=g_vec(gg)
        for sese=1:length(seed_vec)
            seed=seed_vec(sese)
            result=struct;
            rng(seed)
            
            %%measuring drift speed
            dt=0.1;
            tmax_drift=3;
            velocity_sample_timestep=7:8;
            %%net settings
            net.g=0;
            net.N=4000;
            net.phi=@(x)erf(x./sqrt(2));
            net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
            hp=struct;
            sim=struct;
            hp.A=A;
            [hp, net, sim] = prep_network_param(hp, net,sim);
            
            x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
            
            %% learning output weights for classic ring
            sim.r = net.phi(x);
            pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
            regfac = hp.M*hp.alpha_reg*eye(length(pts));
            net.wout = (sim.r(:,pts)/...
                (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
            
            %% obtaining open loop output
            result.baseline.z_ol = net.wout'*sim.r;
            
            
            
            %% adding noise to 'classic' ring
            
            
            net_with_g=net;
            net = update_net_g(net,g);
            x_noi = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
            sim.r = net.phi(x_noi);
            result.with_noise.z_ol = net.wout'*sim.r;
            
            %%
            delta_psi=phase([1,1i]*result.with_noise.z_ol)-sim.psi;
            delta_rho=abs([1,1i]*result.with_noise.z_ol)-hp.A;
            result.rmse=sqrt(mean(delta_psi.^2))
            result.eDC=mean(delta_psi)
            result.rmseNoDc=sqrt(mean(delta_psi.^2)-result.eDC^2)
            %%
            [x_cl,result.with_noise.z_cl]=fast_conv_to_fp_extended(net,[],...
                struct('xinit',x,...
                'ol',0,...
                'tmax',tmax_drift,...
                'dt',dt,'save_neurons',1));
            %%
            result.v=1/dt*mean(diff(mod(angle(result.with_noise.z_cl(:,velocity_sample_timestep)),2*pi),[],2),2);
            %%
            N=net.N;
            save(['HRresultsN',num2str(net.N),'iter',num2str(cnt)],'result','N','g','A','seed');
            cnt=cnt+1;
        end
    end
end