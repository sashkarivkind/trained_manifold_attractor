clear;

result=struct;
tmax=1000;
dt=1;

net.g=1;
hp.omega_vec=0;
%todo: save or force seed for RNG!
[hp, net, sim] = prep_network_param(hp, net,struct);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));


sim.r = net.phi(x);
result.with_learning={};
M_vec=[3,5];
Gij={};
results.G_polar={};
for MM=1:length(M_vec)
    hp.M=M_vec(MM);
    
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
    results.G_polar{MM}=calc_polar_gain(Gij,sim.psi(pts));
end


figure;
for MM=1:2
    G_vec_over_m=[];
    for pp=1:length(results.G_polar{MM}.flat)
        G_vec_over_m(end+1)=results.G_polar{MM}.flat{pp}{2,2};
    end
    plot(G_vec_over_m)
    hold on;
end
