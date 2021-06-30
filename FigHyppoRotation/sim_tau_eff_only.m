   %% simulate

N=4000;
% g_vec=[1];%0:0.1:2;
g_vec=[1];%0:0.1:2;
% h_vec=[0.8 0.2 0 0.2 -0.5];
% h_vec=[-0.15];
h_vec=[-0.2  0.9];
a=4;b=1;
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
                hp.M=40;
                A=1.2;
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
                
                
            save(['TauEffResultsN',num2str(net.N),...
                'a',strrep(num2str(a),'.','p'),...
                'h',strrep(num2str(h),'.','p'),...
                'g',strrep(num2str(net.g),'.','p'),...
                'seed',num2str(seed)]);
            toc;
            end
        end
    end
    
    
end
%% analyze
sim.x=x;
tau_inv_vec=[];
tau_eff2=[];

all_pts=1:length(sim.f_ol);
for ppp=1:length(all_pts)
    semi=semi_empirical_spectrum(net,sim,20,all_pts,ppp);
    % tau_inv_vec(end+1)=semi.tau_inv;
    tau_eff2(end+1)=semi.tau_approx_by_deriv;
    
end

%% visualize tau eff:
figure;
scatter(sim.f_ol(1,:),sim.f_ol(2,:),30,tau_eff2,'o')
colorbar
xlim([-1.5 1.5]);ylim([-1.5 1.5])
set(gcf,'Position',[ 903   775   217   173])
hAx=gca
colormap(hsv)
hAx.CLim=[1 1.5]
lineStyles = linspecer(100,'sequential');
colormap(lineStyles)
hAx.CLim=[0.5 1]