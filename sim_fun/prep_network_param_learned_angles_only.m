function [ hp, net, sim ] = prep_network_param_learned_angles_only(hp_custom,...
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
hp.omega_vec = logspace(-2,4,40);

%setting default network parameters
net.g=1.0;
net.N=1000;
net.phi=@(x) tanh(x); %added 1 to break symmetry
net.phip =@(x) 1-tanh(x).^2;
net.W=net.g*randn(net.N)/sqrt(net.N);


% % % % sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate
% % % % sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)]; This location
% results in BUGS

hp=nom_opt_assigner(hp_custom,hp);
sim.psi=pi*(0:(hp.M-1))/(hp.M); %learned angles to simulate
sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];

net=nom_opt_assigner(net_custom,net);
net = update_net_g(net,net.g);
net.theta=2*pi*(0:net.N-1)/net.N;
net.wfb=[cos(net.theta'),sin(net.theta')];
sim=nom_opt_assigner(sim_custom,sim);
sim.psi=pi*(0:(hp.M-1))/(hp.M); %learned angles to simulate
sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];

end

