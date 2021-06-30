clear;
% net.g=1.0;
% hp.A=0.6
b=0.4;
[hp, net, sim] = prep_network_param(struct,struct,struct);
% net.wfb=randn(net.N,2);
sim.f_ol=diag([1,b])*sim.f_ol;
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
