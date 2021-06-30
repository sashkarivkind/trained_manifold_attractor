close all;clear;
% load('FallRotDt1Yonina_ResultsM_3_N4000a2h0g1eps0p01seed4.mat')
load('FallRotDt1Yonina_ResultsM_3_N16000a2h0g1eps0p01seed1.mat')
x_epsi_on_manifold = fast_conv_to_fp(net,...
    sim.f_ol,...
    struct('ol_with_fixed_input',1,...
    'u',0*u_in(:,1000)*ones(1,size(sim.f_ol,2))));
z_epsi_on_manifold = net.wout'*net.phi(x_epsi_on_manifold);

%%
figure;

for icic=1:n_ic
            hhh=plot(zdeltaFall{icic},'linewidth',5)
            set(hhh,'Color',[get(hhh,'Color'),0.5]);
            hold on
end
xlim([-1.5 1.5]);ylim([-1.5 1.5])
plot(sim.f_ol(1,1:160/6:160),sim.f_ol(2,1:160/6:160),'o','MarkerSize',10)
plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k')

axis square; box off;


figure


figure
psi_on_manifold = atan2(sim.f_ol(2,:),sim.f_ol(1,:));
% plot(psi_on_manifold,atan2(z_epsi_on_manifold(2,:),z_epsi_on_manifold(1,:))-psi_on_manifold )
Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;

plot_with_sorted_x(psi_on_manifold,flipud(Delta))
hold on;
plot_with_sorted_x(psi_on_manifold,zeros(size(psi_on_manifold)),'-')

plot(sim.psi(1:160/6:160)-pi,zeros(6,1),'or','MarkerSize',10)
% axis square; 
box off;xlim([-pi pi]);set(gca,'xtick',[-pi:pi/2:pi])

