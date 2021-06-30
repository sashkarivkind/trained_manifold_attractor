d_spectrum_rr=[];
d_spectrum_pp=[];

for uu=1:length(beta_rescaled_rec)
    v1=eig(M_rescaled_rec{uu});   
    v2=eig(beta_rescaled_rec{uu});

    dist=sort(nn_min_dist(v1,v2),'descend');
    if dd_rec{uu}==1
        d_spectrum_rr(:,end+1)=dist;
    else
        d_spectrum_pp(:,end+1)=dist;
    end
end

figure; 
subplot(2,1,1)
plot(d_spectrum_pp,'x')
title('difference between closed loop and open loop')
ylabel('\Delta\lambda \psi\psi')
xlim([-1,11]);
subplot(2,1,2)
plot(d_spectrum_rr,'x')
ylabel('\Delta\lambda  rr')
xlim([-1,11]);
xlabel('eigne value #')
