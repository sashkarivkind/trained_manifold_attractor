% clear;
grid_res=0.025;
N=16000;
% g_vec=[1];%0:0.1:2;
g_vec=[0.5];%0:0.1:2;
h_vec=[0];
n_ic=10;
DegRot=270;

A=1.2;
a=2;b=1;

clear x_ic
% h_vec=[-0.5 0 0.2 0.4 0.6];
for seed=[1];%8:10;
    seed
    run_new_net = 1;
    zdelta_rec={};
    % g_vec=0:0.2:2;
    input_mode = 'hyppo';
    t_conv_vec=[];
    for gg=1:length(g_vec)
        for hh=1:length(h_vec)
            tic;
            gg
            if 1 %
                rng(seed);
                hp=struct;
                net=struct;
                hp.M=3;
                h=h_vec(hh);
                rot=0;
                sim=struct;
                
                net.g=g_vec(gg);
                net.N=N;
                
                [hp, net, sim] = prep_network_param(hp, net, sim);
                
                sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                
                
                x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
                
                %learning output weights
                sim.r = net.phi(x);
                pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
                regfac = hp.M*hp.alpha_reg*eye(length(pts));
                net.wout = (sim.r(:,pts)/...
                    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
                
                %obtaining open loop output
                sim.z_ol = net.wout'*sim.r;
            end
            net.win=net.wfb;
            
            dt=1;
            tmax=1000;
            epsilonmag=0;
            tswitch=[1];
            pulse_width = [2000]/dt;
            psi_vec=[DegRot]*pi/180;
            nsteps=(round(tmax/dt));
            %spikes
            u_in=zeros(2,nsteps);
            ii=0;
            for ttt =tswitch
                ii=ii+1;
                psi=psi_vec(ii);
                uu1=round(ttt/dt);
                uu2=min(size(u_in,2),uu1+pulse_width(ii));
                %%
                if strcmp(input_mode,'sincos')
                    u_in(:,uu1:uu2)=epsilonmag*[cos(psi);sin(psi)]*...
                        ones(1,size(u_in(:,uu1:uu2),2));
                elseif  strcmp(input_mode,'hyppo')
                    u_in(:,uu1:uu2)=...
                        epsilonmag*1/(a-b)*...
                        [(a-b)*cos(psi-rot)+h*cos((a-b)/b*(psi-rot));...
                        (a-b)*sin(psi-rot)-h*sin((a-b)/b*(psi-rot))]./(1+h/(a-b))...
                        *ones(1,size(u_in(:,uu1:uu2),2));
                else
                    error
                end
                %%
            end
            %             x_ic=2*randn(net.N,n_ic);
           

            x_epsi_on_manifold = fast_conv_to_fp(net,...
                sim.f_ol,...
                struct('ol_with_fixed_input',1,...
                'u',u_in(:,1000)*ones(1,size(sim.f_ol,2))));
            z_epsi_on_manifold = net.wout'*net.phi(x_epsi_on_manifold);
            
            
            save(['FallRotDt1Yonina_ResultsM_3_N',num2str(net.N),...
                'a',strrep(num2str(a),'.','p'),...
                'h',strrep(num2str(h),'.','p'),...
                'g',strrep(num2str(net.g),'.','p'),...
                'eps',strrep(num2str(epsilonmag),'.','p'),...
                'seed',num2str(seed)]);
            toc;
        end
        
        
    end
end