% clear;

% result=struct;
tmax=1000;
dt=1;
g_vec=[1.2:0.2:2];
M_vec=[6];
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;
results.G_polar={};
%
net.N=1000;

for MM=1:length(M_vec)
    hp.M=M_vec(MM);
    MM
    for GG=1:length(g_vec)
        g=g_vec(GG);
        net.g=g;
        GG
        tic;
        %todo: save or force seed for RNG!
        
        for II=1:10
            result.with_learning={}; Gij={};
            [hp, net, sim] = prep_network_param(hp, net,struct);
            hp.sim_resolution=hp.sim_resolution;
            
            x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
            sim.r = net.phi(x);
            pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
            regfac = hp.M*hp.alpha_reg*eye(length(pts));
            net.wout = (sim.r(:,pts)/...
                (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
            
            cnt=0;
            for pt = pts
                cnt=cnt+1;
                [detGxx,Gij{cnt}] = calcGij(net.W*diag(net.phip(x(:,pt)))-eye(net.N),...%internal matrix
                    diag(net.phip(x(:,pt)))*net.wout,...%output
                    net.wfb,... %input
                    hp.omega_vec);
            end
            results.G_polar{MM,GG,II}=calc_polar_gain(Gij,sim.psi(pts));
%             for pp=1:length(results.G_polar{MM,GG,II}.flat)
%                 if(results.G_polar{MM,GG,II}.flat{pp}{2,2}>2)
%                     plot
%                 end
%             end
        end
        toc;
    end
end
%
% save('Gp_10_sim_interations_N_16000_M_12.mat','hp','sim','results')
%%
tmax=1000;
dt=1;
g_vec=[0.2:0.2:2];
% g_vec=[1.2:0.2:2];
% M_vec=[3 5];
g_vec=[1.2 1.4];
M_vec=[6];
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;

figure;
G_vec_over_m_g={};
for MM=1:1
    for GG=1:length(g_vec)
        for II=1:10
            for pp=1:length(results.G_polar{MM,GG,II}.flat)
                G_vec_over_m_g{MM}(GG,pp,II)=results.G_polar{MM,GG,II}.flat{pp}{2,2};
            end
        end
    end
end
%
plot(g_vec,squeeze(G_vec_over_m_g{1}(:,:)),'.r')
hold all
plot(g_vec,mean(G_vec_over_m_g{1}(:,:),2),'-r')
hold on;
% plot(g_vec,squeeze(G_vec_over_m_g{2}(:,:)),'.b')
% plot(g_vec,mean(G_vec_over_m_g{2}(:,:),2),'-b')
