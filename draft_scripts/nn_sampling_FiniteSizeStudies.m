skipsim=1;
if skipsim
clear;

result=struct;
tmax=1000;
dt=1;

net.g=1;
hp.omega_vec=0;
hp.A = 2.0;
%todo: save or force seed for RNG!
[hp, net, sim] = prep_network_param(hp, net,struct);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

end
sim.r = net.phi(x);
result.with_learning={};
M_vec=[5,10,20];
Gij={};
results.G_polar={};
results.G_polar2={};
for MM=1:length(M_vec)
    hp.M=M_vec(MM);
    
    pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
    regfac = hp.M*hp.alpha_reg*eye(length(pts));
    net.wout = (sim.r(:,pts)/...
        (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
    
    CC=sim.r(:,pts)'*sim.r(:,pts);
    idealized_norm_fac = sim.f_ol(1,pts)*CC(:,1);
    net.wout2 = hp.A*sim.r(:,pts)*sim.f_ol(:,pts)'/idealized_norm_fac; %obtain least mean square solution
    
    cnt=0;
    Gij={};
    Gij2={};
    for pt = pts
        cnt=cnt+1;
        [detGxx,Gij{cnt}] = calcGij(net.W*diag(net.phip(x(:,pt)))-eye(net.N),...%internal matrix
            diag(net.phip(x(:,pt)))*net.wout,...%output
            net.wfb,... %input
            hp.omega_vec);
            [detGxx2,Gij2{cnt}] = calcGij(net.W*diag(net.phip(x(:,pt)))-eye(net.N),...%internal matrix
            diag(net.phip(x(:,pt)))*net.wout2,...%output
            net.wfb,... %input
            hp.omega_vec);
    end
    results.G_polar{MM}=calc_polar_gain(Gij,sim.psi(pts));
    results.G_polar2{MM}=calc_polar_gain(Gij2,sim.psi(pts));
end


    figure;
    for MM=1:3
        G_vec_over_m=[];
        for pp=1:length(results.G_polar{MM}.flat)
            G_vec_over_m(end+1)=results.G_polar{MM}.flat{pp}{2,2};
        end
        plot(G_vec_over_m)
        hold on;
    end

    set(gca,'ColorOrderIndex',1);

    for MM=1:3
        G_vec_over_m=[];
        for pp=1:length(results.G_polar2{MM}.flat)
            G_vec_over_m(end+1)=results.G_polar2{MM}.flat{pp}{2,2};
        end
        plot(G_vec_over_m,'--')
        hold on;
    end