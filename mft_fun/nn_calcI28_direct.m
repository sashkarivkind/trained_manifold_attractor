function o=calcI28_direct(crfu,crprp,hmax);
effsize=ceil(hmax/2);
o=zeros(effsize);
psi_vec=2*pi*(0:(length(crprp)-1))/length(crprp);
for nn=1:effsize
    for mm=1:effsize
        nnt=nn*2-1;%taking odd coefficients
        mmt=mm*2-1;
        o(nn,mm)=...
            crfu(mm*2)*...
            2*mean(crprp.*sin(nnt*psi_vec).*sin(mmt*psi_vec));
    end
end