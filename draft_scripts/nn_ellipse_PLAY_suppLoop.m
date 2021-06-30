saves={'ellipse_g1p0_M40_h2side_N2K_fixedNwk',
    'ellipse_g1p0_M40_h2side_N2K_fixedNwk2',
    'ellipse_g1p0_M40_h2side_N2K_fixedNwk3'
    };
for ff=1:length(saves)
    load(saves{ff});
for hh=1:length(h_vec)

    for kk=1:length(K_vec)
        K=K_vec(kk);
        Grr(hh,kk)=o{hh,kk}.G(1,1); 
        Gpp(hh,kk)=o{hh,kk}.G(2,2); 
        wcos(:,hh)=wo_rec{hh}(:,1);
% EErad=o{hh,kk}.G(2,2)
    end
            thiso=o{hh,2};
        MM1=inv(thiso.lhs)* (thiso.B1+thiso.B0(:,1)*((thiso.qn(1,:)*diag(thiso.Gnprefac))));
        eemax(hh)=max(real(eig(MM1)));
end
figure(139); %plot(h_vec,Grr(:,end),'b');hold on;
% figure(38); plot(h_vec,sqrt(mean(wcos.^2,1))); hold on;
 plot(h_vec,eemax,'m'); 
end
% %%
% kkk=(1:10)';
% % U=[cos(kkk*sim.psi(pts));sin(kkk*sim.psi(pts))]*o{hh,kk}.eta;
% U=[cos(kkk*sim.psi(pts));sin(kkk*sim.psi(pts))]*zzz.eta;