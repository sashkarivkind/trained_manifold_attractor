clear;
tmax=1000;

A_vec=[0.5,1,1.2,1.5,2.];
J2_vec=[];
cnt=0;
for aa=1:length(A_vec)
    A=A_vec(aa);
    result=struct;
    
    %%measuring drift speed
    dt=0.1;
    tmax_drift=3;
    velocity_sample_timestep=7:8;
    %%net settings
    net.g=0;
    net.N=1000;
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
    J2_vec(end+1)=trace(net.wout'*net.wfb)
end
