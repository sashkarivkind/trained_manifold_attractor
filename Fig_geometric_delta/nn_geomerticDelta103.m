clear;

epsilon=1e-2;
psi_vec=pi*(0:0.2:0.9999);

net.N=500;
net.g=1.0;
hp.A=1.2
[hp, net, sim] = prep_network_param(hp,net,struct);


a=2;b=1; rot=0;
h=0.6;
sim.f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
psi_emp=phase2(sim.f_ol);


x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
sim.x=x;
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
%% calculating open loop gain
vv={};
GG={};
phaseG=[];
for ptpt = 1:length(sim.psi)
    semi=semi_empirical_spectrum(net,sim,20,1:length(sim.psi),ptpt);
    [vv{ptpt},GG{ptpt}]=eig(semi.G);
    [~,ii] = sort(diag(GG{ptpt}));
    vv{ptpt} = vv{ptpt}(:,ii);
    phaseG(ptpt,:)=phase2(vv{ptpt});
 end

%%
%%
dz=diff(sim.f_ol(:,[end,1:end]),1,2); %differntial of z 
dz_norm = sqrt(sum(dz.^2));
% dpsi=diff(sim.psi([end,1:end])); %differntial of z 
% dpsi=diff(phase2(sim.f_ol(:,[end,1:end])),1,2); %differntial of z 
dpsi=diff(psi_emp(:,[end,1:end]),1,2); %differntial of z 
dpsi(1)=dpsi(2); %compenstaing for discontinuity
dz_dpsi = dz_norm./dpsi;
emp_tan=phase2(diff(sim.f_ol(:,[end,1:end]),1,2)); %tangential angle
emp_norm=emp_tan-pi/2;

%%
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
    delta_z_geom{psi_psi}=f_geom{psi_psi}-sim.f_ol;
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
    plot(psi_emp,delta_epsi_full{psi_psi},'x')
    set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
    plot(psi_emp,delta_epsi_geom{psi_psi},'o')
    set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);

    
    for uu=1:length(sim.psi)
        delta_z_in_ev = inv(vv{uu})*delta_z_geom{psi_psi}(:,uu)
        delta_z_geom_par{psi_psi}(uu)=delta_z_in_ev(2)/dz_dpsi(uu);
%         delta_z_geom_perp{psi_psi}=f_geom{psi_psi}-sim.f_ol;
    end
    plot(psi_emp,delta_z_geom_par{psi_psi},'*')


end


%%
figure; 
plot(180/pi*psi_emp, psi_emp - emp_norm);
hold on;
plot(180/pi*psi_emp(2:end),mod(phaseG(2:end,1),pi)-mod(psi_emp(2:end),pi)');

% plot(sim.psi, sim.psi - emp_norm);
%%
figure; 
plot(psi_emp, dz_dpsi);
%%
%     figure;
%     plot(phase2(sim.z_ol)-phase2(sim.f_ol));

%ensuring that phase of G-parallel as obtained from e.v. is same
%as tangential direction from manifold geometry... 
%checking that the differnce is zero up to finite resolution on psi
% figure;
% plot(mod(phaseG(2:end,2),pi)-mod(emp_tan(2:end),pi)');
% hold on;
% plot(mod(phaseG(1:end-1,2),pi)-mod(emp_tan(2:end),pi)');
% %%
% figure;
% plot(mod(phaseG(:,1),pi)-mod(emp_norm,pi)');
% % hold on;
% % plot(mod(phaseG(:,1),pi)-mod(emp_norm,pi)');
% 
% %%
% figure;
% % plot(mod(phaseG(:,1),pi)-mod(sim.psi,pi)');
% % hold on;
% % plot(mod(phaseG(:,1),pi)-mod(emp_norm,pi)');
% % plot(mod(phaseG(:,1),pi)-mod(psi_emp,pi)');
% plot(180/pi*(mod(phaseG(:,1),pi)-mod(psi_emp,pi)'));
% hold on;
% plot(180/pi*(mod(phaseG(:,1),pi)-mod(emp_norm,pi)'));
% 
% plot(180/pi*(mod(phaseG(:,2),pi)-mod(emp_tan,pi)'));
% 
% % hold on;
% % plot(mod(phaseG(1:end-1,1),pi)-mod(emp_norm,pi)');