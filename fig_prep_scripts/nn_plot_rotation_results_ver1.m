function nn_plot_rotation_results_ver1(sim,hp,result,psi,epsilonmag,title_pfx)
if nargin<6
   title_pfx='';
end
figure(100);
xaxis_deg=(0:(size(sim.f_ol,2)-1))*2*pi/size(sim.f_ol,2);
subplot(2,2,1);
plot(xaxis_deg,phase([1,1i]*result.rotation.zdelta_ol)-sim.psi)
title([title_pfx,'rotation']);
figure(101);
plot(angle(result.rotation.zdelta))
ylabel('transient');

