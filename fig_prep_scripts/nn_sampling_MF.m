function [G,debuu]=nn_sampling_MF(phi,J_ext,J_0,Amp,K,psi_vec,phip,omega_vec)
debuu={};
[mu_out, sig_out,beta_out,beta_add]=...
    solveFixedPointRing_Erf(phi,J_ext,J_0,Amp,1,1,K,psi_vec);
sig_out=sqrt(sig_out.^2-Amp^2/2);
beta_out=beta_out-Amp.^2.*beta_add;

crr=crr_fun(sig_out,beta_out,psi_vec,Amp,phi);

cc=beta_out/sig_out.^2;
hmax=2;
cfu=fft(crr);
cfu=2*real(cfu)/length(cfu); %cosine fourier coefficients
