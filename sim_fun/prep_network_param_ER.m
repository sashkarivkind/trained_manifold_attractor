function [ hp, net, sim ] = prep_network_param_ER(hp_custom,...
    net_custom,...
    sim_custom)
%PREP_NETWORK_PARAM Summary of this function goes here
%   Detailed explanation goes here

%setting default hyper parameters
hp.alpha_reg = 0.0;
hp.b=1;
hp.M=20;
hp.A=1.2;
hp.sim_resolution=160; %points to simulate on the ring

net=net_custom; %% UPDATE THIS PART!!!
%setting default network parameters
% net.K=400;
% net.g=1.0;
% net.N=1000;
p=net.K/net.N;
net.W  = p>rand(net.N,net.N);
net.W=net.W-diag(diag(net.W));
net.W =sparse(net.W );
% andd somewhere INP sqrt(K)

net.phi=@(x) (x).*(x>0);
net.phip =@(x) 1.*(x>0);


% % % % sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate
% % % % sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)]; This location
% results in BUGS

hp=nom_opt_assigner(hp_custom,hp);
sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate WHY ITS DOUBLE??
sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];

net=nom_opt_assigner(net_custom,net);

net.W = net.W.*(-net.g/sqrt(net.K));

% net = update_net_g(net,net.g);
net.theta=2*pi*(0:net.N-1)/net.N;
net.wfb=[cos(net.theta'),sin(net.theta')];
sim=nom_opt_assigner(sim_custom,sim);
sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate
sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];

end

