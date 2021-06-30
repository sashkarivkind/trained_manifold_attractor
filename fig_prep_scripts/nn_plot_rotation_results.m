function nn_plot_rotation_results(sim,hp,result,psi,epsilonmag,title_pfx)
if nargin<6
   title_pfx='';
end
figure;
xaxis_deg=(0:(size(sim.f_ol,2)-1))*360/size(sim.f_ol,2);
subplot(2,2,1);
plot(xaxis_deg,phase([1,1i]*result.rotation.zdelta_ol)-sim.psi)
title([title_pfx,'rotation']);

hold on; 
plot(xaxis_deg,phase([1,1i]*result.rotation.zdelta_ol_ref)-sim.psi)
grid
ylabel('\Delta\psi');
subplot(2,2,2);
plot(xaxis_deg,abs([1,1i]*result.rotation.zdelta_ol)-hp.A)
hold on; 
plot(xaxis_deg,abs([1,1i]*result.rotation.zdelta_ol_ref)-hp.A)
grid
ylabel('\Delta\rho');

subplot(2,2,3);
plot(angle(result.rotation.zdelta)/pi*180)
ylabel('transient');


subplot(2,2,4);
fol=[1,1i]*sim.f_ol;
fol_new= fol+epsilonmag*(cos(psi)+1i*sin(psi));
xaxis_deg=(0:(size(sim.f_ol,2)-1))*360/size(sim.f_ol,2);

plot(xaxis_deg,phase([1,1i]*result.rotation.zdelta_ol)-phase([1,1i]*result.rotation.zdelta_ol_ref))
hold on;
plot(xaxis_deg,phase(fol_new)-sim.psi)
ylabel('error in/out');
