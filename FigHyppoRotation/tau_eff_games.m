
sim.x=x;
tau_inv_vec=[];
tau_eff2=[];

all_pts=1:length(sim.f_ol);
for ppp=1:length(all_pts)
semi=semi_empirical_spectrum(net,sim,20,all_pts,ppp);
% tau_inv_vec(end+1)=semi.tau_inv;
tau_eff2(end+1)=semi.tau_approx_by_deriv;

end

%%
%%
figure
st=20;
for icic=1:n_ic
    psi_dot=diff(angle(zdeltaRot{icic}))/dt;

            plot(angle(zdeltaRot{icic}(st-1+2:end)),0.0000*icic+psi_dot(st:end),'.')
            hold on
end
% plot(psi_on_manifold,tau_inv_vec.*Delta,'k-')
psi_on_manifold = atan2(sim.f_ol(2,:),sim.f_ol(1,:));
% plot(psi_on_manifold,atan2(z_epsi_on_manifold(2,:),z_epsi_on_manifold(1,:))-psi_on_manifold )
Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;

plot(psi_on_manifold,Delta./tau_eff2,'k-')
plot(psi_on_manifold,Delta,'k:')
% plot(psi_on_manifold,Delta,'k--')
axis square; box off;
% xlim([0 1000]);

xlim([-pi pi]);ylim([-1 1]*0.01);set(gca,'Xtick',[-pi:pi/4:pi])

%%
figure;
% plot(psi_on_manifold,tau_inv_vec.*Delta)
plot(psi_on_manifold,Delta./tau_eff2)
hold on;
plot(psi_on_manifold,Delta,'linewidth',2)


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
% %  contour(z_vec,z_vec,reshape(dpsiu,size(z_grid)),-1:1,'r--');
contour(z_vec,z_vec,reshape(dpsiu,size(z_grid)),-5e-2:5e-3:5e-2,'r--');
hold on
contour(z_vec,z_vec,reshape(rhou-rho_aux,size(z_grid)),-1:1:1,'b--');
plot(z_epsi_on_manifold(1,:),z_epsi_on_manifold(2,:),'k')
% plot(real(zdelta(1:10:end)),imag(zdelta(1:10:end)),'og','MarkerSize',5);

%%
figure;
scatter(sim.f_ol(1,:),sim.f_ol(2,:),60,tau_eff2,'*')
box off;axis square
xlim([-1.5 1.5]);ylim([-1.5 1.5]);
colorbar
