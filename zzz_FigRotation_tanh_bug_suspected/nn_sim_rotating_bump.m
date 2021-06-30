run_new_net = 1;
if run_new_net
clear;
rng(1);
hp=struct;
net=struct;
sim = struct;
net.g=1.6;
[hp, net, sim] = prep_network_param(hp, net, sim);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%learning output weights
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

%simulating closed loop
x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
r_test_rand=net.phi(x_test_rand);
sim.z_test_rand= net.wout'*r_test_rand;

x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
r_test_noi=net.phi(x_test_noi);
sim.z_test_noi= net.wout'*r_test_noi;

%plotting results
figure;
plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
hold on;
plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);

end
net.win=net.wfb;

dt=0.1;
tmax=150;
epsilonmag=0.3;
tswitch=[10,40,80];
psi_vec=[95,30,45]*pi/180;
nsteps=(round(tmax/dt));
%spikes
u_in=zeros(2,nsteps);
ii=0;
for ttt =tswitch
    ii=ii+1;
    psi=psi_vec(ii);
u_in(:,round(ttt/dt):end)=epsilonmag*[cos(psi);sin(psi)]*...
    ones(1,size(u_in(:,round(ttt/dt):end),2));
end
%tracking
% dpsi_dt=0.1;
% t_vec=dt:dt:tmax;
% u_in=[-sin(t_vec*dpsi_dt) ;cos(t_vec*dpsi_dt) ];
[xdelta,zdelta,xrecdelta]=fast_conv_to_fp_extended(net,[],...
    struct('xinit',x(:,1),...
    'ol',0,...
    'tmax',tmax,...
    'dt',dt,...
    'u_in',u_in,...
    'save_neurons',[1:50:1000]));

figure; 
subplot(2,1,1);
plot(angle(zdelta)); 
hold on;
plot(abs(zdelta))
subplot(2,1,2);
plot(xrecdelta');
