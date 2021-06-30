clear;

epsilon=1e-2;
psi_vec=pi*(0:0.2:0.9999);

net.N=500;
net.g=1.0;
hp.A=1.5
[hp, net, sim] = prep_network_param(hp,net,struct);


a=2;b=1; rot=0;
h=0.5;
sim.f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));



x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

u={};
z_epsi={};
for psi_psi = 1:length(psi_vec)
    net.win=net.wfb;
    u{psi_psi}=epsilon*diag([cos(psi_vec(psi_psi));sin(psi_vec(psi_psi))])...
        *ones(size(sim.f_ol));
    f_ol_tag = sim.f_ol + u{psi_psi};
    x_epsi = fast_conv_to_fp(net,...
        sim.f_ol,...
        struct('ol_with_fixed_input',1,...
        'u',u{psi_psi}));
    z_epsi{psi_psi}=net.wout'*net.phi(x_epsi);
    f_geom{psi_psi}=sim.f_ol+u{psi_psi};
    delta_epsi_full{psi_psi}=phase2(z_epsi{psi_psi})-phase2(sim.z_ol);
    delta_epsi_geom{psi_psi}=phase2(f_geom{psi_psi})-phase2(sim.f_ol);
    psi_psi
end

x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
r_test_noi=net.phi(x_test_noi);
sim.z_test_noi= net.wout'*r_test_noi;

%% plotting results
figure;
subplot(1,2,1)
plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
hold on;
plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);



subplot(1,2,2)
hold on;
for psi_psi = 1:length(psi_vec)
    plot(sim.psi,delta_epsi_full{psi_psi},'x')
    set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
    plot(sim.psi,delta_epsi_geom{psi_psi},'o')
end

%%
%     figure;
%     plot(phase2(sim.z_ol)-phase2(sim.f_ol));