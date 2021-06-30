clear;
result={};
g_vec=[0.01 0.5,1, 1.5,1.8];
for gg=1:length(g_vec)
    
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
net.g=g_vec(gg);
[hp, net, sim] = prep_network_param(struct,net,struct);

x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%learning output weights
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

result{gg}.C=(sim.r'*sim.r)/net.N/hp.sim_resolution;
end
%plotting results
figure;
% subplot(2,1,1)
% plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
% hold on;
% plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
% subplot(2,1,2)
for gg=1:5
semilogy(flip(eig(result{gg}.C)));
hold on;
end