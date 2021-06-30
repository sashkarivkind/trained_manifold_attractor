function [detGxx,Gij] = calcGij(W,vv,uu,omega_vec)
    %Wprime=W+uu*vv'; so uu is feedback and vv is the readout
    for kk=1:length(omega_vec)
        Gij{kk} = vv'/(1i*omega_vec(kk)*eye(size(W))-W)*uu;
        detGxx(kk) = det(Gij{kk}-eye(size(Gij{kk})));
    end