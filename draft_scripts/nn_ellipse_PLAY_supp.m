for hh=1:length(h_vec)

    for kk=1:length(K_vec)
        K=K_vec(kk);
        Grr(hh,kk)=o{hh,kk}.G(1,1); 
        Gpp(hh,kk)=o{hh,kk}.G(2,2); 
% EErad=o{hh,kk}.G(2,2)
        thiso=o{hh,kk};
        MM1=inv(thiso.lhs)* (thiso.B1+thiso.B0(:,1)*((thiso.qn(1,:)*diag(thiso.Gnprefac))));
        eemax(hh,kk)=max(real(eig(MM1)));
    end
        wcos(:,hh)=wo_rec{hh}(:,1);

end
figure; plot(h_vec,Grr(:,end))
figure; plot(h_vec,eemax(:,end))
% %%
% kkk=(1:10)';
% % U=[cos(kkk*sim.psi(pts));sin(kkk*sim.psi(pts))]*o{hh,kk}.eta;
% U=[cos(kkk*sim.psi(pts));sin(kkk*sim.psi(pts))]*zzz.eta;