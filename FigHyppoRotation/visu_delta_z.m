% clear z_epsi0 z_epsi
% z_vec=0:0.025:2;
% rr=1;
% z_grid = repmat(z_vec,length(z_vec),1);
% z1=z_grid(:);
% z_grid_t=z_grid';
% z2=z_grid_t(:);
% f_ol_tag=[z1';z2'];
% psi_aux=atan2(f_ol_tag(2,:),f_ol_tag(1,:));
% rho_aux=abs(f_ol_tag(2,:)+1i*f_ol_tag(1,:));
% 
% x_epsi0 = fast_conv_to_fp(net,...
%     f_ol_tag,...
%     struct('ol_with_fixed_input',1,...
%     'u',0*u_in(:,1000)*ones(1,size(f_ol_tag,2))));
% z_epsi0{rr} = net.wout'*net.phi(x_epsi0);
% 
% x_epsi = fast_conv_to_fp(net,...
%     f_ol_tag,...
%     struct('ol_with_fixed_input',1,...
%     'u',u_in(:,1000)*ones(1,size(f_ol_tag,2))));
% z_epsi{rr} = net.wout'*net.phi(x_epsi);
% %%
% x_epsi0p1 = fast_conv_to_fp(net,...
%     f_ol_tag,...
%     struct('ol_with_fixed_input',1,...
%     'u',0.1*u_in(:,1000)*ones(1,size(f_ol_tag,2))));
% z_epsi0p1{rr} = net.wout'*net.phi(x_epsi0p1);
% 

x_epsi_on_manifold = fast_conv_to_fp(net,...
    sim.f_ol,...
    struct('ol_with_fixed_input',1,...
    'u',u_in(:,1000)*ones(1,size(sim.f_ol,2))));
z_epsi_on_manifold = net.wout'*net.phi(x_epsi_on_manifold);


%%
d0=z_epsi0-f_ol_tag;
d0abs=abs([1,1i]*d0);
psi_0=angle([1,1i]*z_epsi0);
rho_0=abs([1,1i]*z_epsi0);

% d0p1=z_epsi0p1-f_ol_tag;
% d0p1abs=abs([1,1i]*d0p1);
% psi_0p1=angle([1,1i]*z_epsi0p1);
% rho_0p1=abs([1,1i]*z_epsi0p1);

du=z_epsi-f_ol_tag;
duabs=abs([1,1i]*du);
psiu=angle([1,1i]*z_epsi);
rhou=abs([1,1i]*z_epsi);

ddu = du-d0;

%%
% figure;
% quiver(f_ol_tag(1,:),f_ol_tag(2,:),d0(1,:),d0(2,:));
% %%
% figure;
% quiver(f_ol_tag(1,:),f_ol_tag(2,:),du(1,:),du(2,:));
% %%
% %%
% figure;
% quiver(f_ol_tag(1,:),f_ol_tag(2,:),ddu(1,:),ddu(2,:));
%%
figure;
d0abs_zz=d0abs;
d0abs_zz(1)=d0abs_zz(2);
imagesc(z_vec,z_vec,reshape(log10(d0abs_zz),size(z_grid)),[-3,0]);%%
set(gca,'YDir','normal')
xlim([-1.5 1.5]);ylim([-1.5 1.5])

%%
figure;
d0abs_zz=d0abs;
d0abs_zz(1)=d0abs_zz(2);
imagesc(z_vec,z_vec,reshape(log10(d0abs_zz),size(z_grid)),[-3,0]);%%
set(gca,'YDir','normal')
hold on
for icic=1:n_ic
            plot(zdeltaFall{icic},'o')
            hold on
end
xlim([-1.5 1.5]);ylim([-1.5 1.5])
axis square; box off;
%%
figure
for icic=1:n_ic
            plot([1:length(angle(zdeltaRot{icic}))].*dt,angle(zdeltaFall{icic}),'.')
            hold on
end
axis square; box off;
xlim([0 1000]);ylim([-pi pi])
%%
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% figure
% hold on
% plot(hp.A*cos(sim.psi(1:hp.sim_resolution/(2*hp.M):end)),hp.A*sin(sim.psi(1:hp.sim_resolution/(2*hp.M):end)),'or','MarkerSize',10)
% axis square

%%
DeltaAng=0.56;
figure
psi_on_manifold = atan2(sim.f_ol(2,:),sim.f_ol(1,:));
% plot(psi_on_manifold,atan2(z_epsi_on_manifold(2,:),z_epsi_on_manifold(1,:))-psi_on_manifold )
Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;

plot(psi_on_manifold,Delta)
hold on;
plot(psi_on_manifold,zeros(size(psi_on_manifold)),'-')
plot(DeltaAng,0,'o','MarkerSize',10)
axis square; box off;
%%

figure;
imagesc(z_vec,z_vec,reshape(log10(duabs),size(z_grid)),[-3,0]);
set(gca,'YDir','normal')
hold on
axis square; box off;
xlim([-1.5 1.5]);ylim([-1.5 1.5])
colormap('gray')
%%
figure
pointsize = 20;
scatter(sim.f_ol(1,:),sim.f_ol(2,:), pointsize, log10(abs(Delta)));
% set(gca,'clim',[-3 0])
%%
figure;
imagesc(z_vec,z_vec,reshape(log10(duabs),size(z_grid)),[-3,0]);
set(gca,'YDir','normal')
colormap('gray')

hold on
for icic=1:n_ic
            plot(zdeltaRot{icic},'o','markersize',2)
            hold on
end
xlim([-1.5 1.5]);ylim([-1.5 1.5])
axis square; box off;
ang=find(abs(sim.psi-DeltaAng)<(2*pi/160));
plot(sim.f_ol(1,ang(1)),sim.f_ol(2,ang(1)),'or','MarkerSize',10)
%%
figure
for icic=1:n_ic
            plot([1:length(angle(zdeltaRot{icic}))].*dt,angle(zdeltaRot{icic}),'.')
            hold on
end
axis square; box off;
xlim([0 1000]);ylim([-pi pi])

%%
figure
for icic=1:n_ic
    psi_dot=diff(angle(zdeltaRot{icic}))/dt;

            plot(angle(zdeltaRot{icic}(2:end)),0.0000*icic+psi_dot,'.')
            hold on
end
plot(psi_on_manifold,Delta,'k-')
axis square; box off;
% xlim([0 1000]);
ylim([-1 1]*0.01)

%%
% figure;
% imagesc(reshape(log10(d0p1abs),size(z_grid)));
% set(gca,'YDir','normal')

%%
% % figure;
% % imagesc(reshape(ddu(1,:),size(z_grid)));
% % set(gca,'YDir','normal')
% % %%
% % figure;
% % imagesc(reshape(psi_aux,size(z_grid)));
% % set(gca,'YDir','normal')

%%
% figure;
% imagesc(reshape(psi_0p1,size(z_grid)));
% set(gca,'YDir','normal')
%%
% figure;
% dpsi0p1=psi_0p1-psi_aux;
% imagesc(reshape(dpsi0p1,size(z_grid)));
% set(gca,'YDir','normal','CLim',[-0.003,0.003])
%%
% figure;
% imagesc(reshape(rho_0p1-rho_aux,size(z_grid)));
% set(gca,'YDir','normal','CLim',[-0.1,0.1])
%%
% figure;
% dpsi0p1=psi_0p1-psi_aux;
% contour(reshape(dpsi0p1,size(z_grid)),-0.004:0.002:0.004);
% hold on
% contour(reshape(rho_0p1-rho_aux,size(z_grid)),-0.02:0.01:0.02);
%%
figure;
% dpsi0p1=psi_0p1-psi_aux;
dpsi0=psi_0-psi_aux;
dpsiu=psiu-psi_aux;
% contour(z_vec,z_vec,reshape(dpsi0p1,size(z_grid)),-1:1:1,'r');
% hold on
% contour(z_vec,z_vec,reshape(rho_0p1-rho_aux,size(z_grid)),-1:1:1,'b');

% figure;
contour(z_vec,z_vec,reshape(dpsi0,size(z_grid)),-1:1:1,'r:');
hold on
contour(z_vec,z_vec,reshape(rho_0-rho_aux,size(z_grid)),-1:1:1,'b:');
% figure;
contour(z_vec,z_vec,reshape(dpsiu,size(z_grid)),-1:1:1,'r--');
hold on
contour(z_vec,z_vec,reshape(rhou-rho_aux,size(z_grid)),-1:1:1,'b--');
plot(z_epsi_on_manifold(1,:),z_epsi_on_manifold(2,:),'k')
plot(real(zdelta(1:10:end)),imag(zdelta(1:10:end)),'og','MarkerSize',5);

figure
psi_on_manifold = atan2(sim.f_ol(2,:),sim.f_ol(1,:));
plot(psi_on_manifold,atan2(z_epsi_on_manifold(2,:),z_epsi_on_manifold(1,:))-psi_on_manifold )

% set(gca,'YDir','normal','CLim',[-0.003,0.003])
% figure;
% 
% for rr=1:length(r_vec)
% plot(180/pi*sim.psi,angle(z_epsi0(1,:)+1i*z_epsi0(2,:))-sim.psi,'x-')
% hold on
% set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
% plot(180/pi*sim.psi,angle(z_epsi(1,:)+1i*z_epsi(2,:))-sim.psi,'o-')
% rr
% end
% grid
%%
% % figure;
% % plot(real(zdelta),imag(zdelta),'o');
% figure
% for icic=1:n_ic
%             plot(zdeltaRot{icic})
%             hold on
% end

%%
f_ol_upd=sim.f_ol+u_in(:,1000)*ones(1,size(sim.f_ol,2));
Delta_geom = atan2(f_ol_upd(2,:),f_ol_upd(1,:)) - psi_on_manifold;
figure;
plot(Delta,'x');
hold on;
plot(Delta_geom,'o')