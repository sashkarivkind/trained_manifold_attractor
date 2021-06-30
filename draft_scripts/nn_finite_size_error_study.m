clear;
hp.M = 40;
if 1
[hp, net, sim] = prep_network_param(hp,struct,struct);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%learning output weights
sim.r = net.phi(x);
end
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));

net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

%simulating closed loop
% x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
% r_test_rand=net.phi(x_test_rand);
% sim.z_test_rand= net.wout'*r_test_rand;
% 
% x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
% r_test_noi=net.phi(x_test_noi);
% sim.z_test_noi= net.wout'*r_test_noi;

% %plotting results
% figure;
% plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
% hold on;
% plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
% plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
% plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
% plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
%%
psi_ol=phase([1 1i]*sim.z_ol);
rho_ol=abs([1 1i]*sim.z_ol);
C=1/net.N/hp.sim_resolution*sim.r'*sim.r;
[vv,ee]=eig(C);
[ee,ii]=sort(real(diag(ee)),'descend');
vv=real(vv(:,ii));
vvrr=sim.r*vv*diag(real(sqrt(1./ee)));
wowo=vvrr'*net.wout;
%%
figure;
plot(sim.psi,psi_ol-sim.psi)
hold on;
plot(sim.psi,rho_ol-hp.A)
%%
figure;
semilogy(abs(wowo),'x');
%%
figure;
k=10;
for kk=1:k
    subplot(k,1,kk);
    plot(vv(:,2*kk-1));
    hold on
    plot(vv(:,2*kk));
end