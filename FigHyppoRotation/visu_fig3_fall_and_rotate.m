clear;
load('Hypp30DegresultsN1000h-0p15g1seed1_Fall_Rot.mat')

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
colormap('gray');
colorbar;
xlim([-1.5 1.5]);ylim([-1.5 1.5])

%%
figure;
d0abs_zz=d0abs;
d0abs_zz(1)=d0abs_zz(2);
imagesc(z_vec,z_vec,reshape(log10(d0abs_zz),size(z_grid)),[-3,0]);%%
set(gca,'YDir','normal')
colormap('gray');
colorbar;
hold on
for icic=1:n_ic
            hhh=plot(zdeltaFall{icic},'linewidth',5)
            set(hhh,'Color',[get(hhh,'Color'),0.5]);
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
plot(DeltaAng,0,'o','MarkerSize',10,'linewidth',3)
axis square; box off;
%%

figure;
imagesc(z_vec,z_vec,reshape(log10(duabs),size(z_grid)),[-3,0]);
set(gca,'YDir','normal')
colormap('gray');
colorbar;

hold on
axis square; box off;
xlim([-1.5 1.5]);ylim([-1.5 1.5])

%%
figure;
imagesc(z_vec,z_vec,reshape(log10(duabs),size(z_grid)),[-3,0]);
set(gca,'YDir','normal')
colormap('gray')
colorbar;
hold on
for icic=1:n_ic
%             plot(zdeltaRot{icic},'o','markersize',2)
            hhh=plot(zdeltaRot{icic},'linewidth',5)
            set(hhh,'Color',[get(hhh,'Color'),0.5]);
            hold on
end
xlim([-1.5 1.5]);ylim([-1.5 1.5])
axis square; box off;
ang=find(abs(sim.psi-DeltaAng)<(2*pi/160));
% plot(sim.f_ol(1,ang(1)),sim.f_ol(2,ang(1)),'or','MarkerSize',10,'linewidth',3)
plot(sim.f_ol(1,ang(1)),sim.f_ol(2,ang(1)),'or','MarkerSize',10)
%%
figure
for icic=1:n_ic
            plot([1:length(angle(zdeltaRot{icic}))].*dt,angle(zdeltaRot{icic}),'.')
            hold on
end
axis square; box off;
xlim([0 1000]);ylim([-pi pi])


