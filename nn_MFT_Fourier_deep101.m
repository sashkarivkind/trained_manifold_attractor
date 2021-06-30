function [G,debuu,sfpro]=nn_MFT_Fourier_deep101(phi,...
    J_ext,...
    J_0,...
    Amp,...
    K,...
    psi_vec,...
    phip,...
    omega_vec,...
    sfpr,...
    hmax...
)
debuu={};

if nargin < 9
    dosfpr=1;
elseif sfpr==0
    dosfpr=1;
else
    dosfpr=0;
end

if nargin<10
    hmax=8;
end

if dosfpr
    sfpr=struct; %sfpr stands for solveFixedPointRing
    [sfpr.mu_out, sfpr.sig_out,sfpr.beta_out,sfpr.beta_add]=...
        solveFixedPointRing_Erf(phi,J_ext,J_0,Amp,1,1,K,psi_vec);
end

sfpro=sfpr; %saving a copy here, before manipulations on various fields
sfpr.sig_out=sqrt(sfpr.sig_out.^2-Amp^2/2);
sfpr.beta_out=sfpr.beta_out-Amp.^2.*sfpr.beta_add;

crr=crr_fun(sfpr.sig_out,sfpr.beta_out,psi_vec,Amp,phi);
crprp=crr_fun(sfpr.sig_out,sfpr.beta_out,psi_vec,Amp,phip);

cfu=fft(crr);
cfu=2*real(cfu)/length(cfu); %cosine fourier coefficients

cpfu=fft(crprp);
cpfu=2*real(cpfu)/length(cpfu); %cosine fourier coefficients

cfu_direct=crr_fourier_fun(sigma,g,A,phi);


[BETA0,BETA1]=nn_phiphpfourierMC_even(phi,phip,J_0^2*cfu(2:1:hmax),Amp); %todo - align with cfu&CxX
[BETA0pp,BETA1pp]=nn_phiphpfourierMC_evenPsiPsi(phi,phip,J_0^2*cfu(2:1:hmax),Amp); %todo - align with cfu&CxX

% BETA1ppI28=2*J_0^2*calcI28(cfu,cpfu,hmax);
cxfu=J_0^2*cfu; %correlation function for x
BETA1ppI28=nn_calcI28_direct(cxfu,crprp,hmax);
CxX=J_0^2*cfu(2:2:hmax);

Q=J_0^2*inv(diag(CxX))*BETA1pp;
Amp*J_0^2*BETA0pp(1,1)./(1-Q(1,1))./CxX(1)% 1=1

G.rr=zeros(size(omega_vec));
G.pp=zeros(size(omega_vec));
G.rr_approx=zeros(size(omega_vec));
G.pp_approx=zeros(size(omega_vec));
for ii = 1:length(omega_vec)
    iw=1i*omega_vec(ii);
    alpha_vec=inv((1+iw)*diag(CxX)-J_0^2*BETA1)*J_0^2*BETA0;
    alphapp_vec=inv((1+iw)*diag(CxX)-J_0^2*BETA1pp)*J_0^2*BETA0pp;
    G.rr(ii)=Amp*alpha_vec(1);
    G.pp(ii)=Amp*alphapp_vec(1);
    
    Qp=J_0^2*inv(diag(CxX))*BETA1pp;
    G.pp_approx(ii)=Amp*J_0^2*BETA0pp(1,1)./(1+iw-Qp(1,1))./CxX(1);
    Qr=J_0^2*inv(diag(CxX))*BETA1;
    G.rr_approx(ii)=Amp*J_0^2*BETA0(1,1)./(1+iw-Qr(1,1))./CxX(1);
end

debuu.pp_p1.v1=Amp*J_0^2*BETA0pp/CxX(1);
debuu.pp_p1.v2=1-J_0^2*BETA1pp/CxX(1);
debuu.rr_p1.v1=Amp*J_0^2*BETA0/CxX(1);
debuu.rr_p1.v2=1-J_0^2*BETA1/CxX(1);
debuu.betarr={BETA0,BETA1};
debuu.betapp={BETA0pp,BETA1pp};
debuu.betapp128=BETA1ppI28;

beta_closedloop=@(B0,B1) inv(diag(CxX))*(J_0^2*B0*Amp...
    *[1:length(B0)==1]+J_0^2*B1-diag(CxX));
debuu.spectrum_pp=beta_closedloop(BETA0pp,BETA1pp);
debuu.spectrum_rr=beta_closedloop(BETA0,BETA1);
debuu.crr=crr;
debuu.crprp=crprp;

%%

for ii = 1:length(omega_vec)
    iw=1i*omega_vec(ii);
    alpha_vec=inv((1+iw)*diag(CxX)-J_0^2*BETA1)*J_0^2*BETA0;
    alphapp_vec=inv((1+iw)*diag(CxX)-J_0^2*BETA1pp)*J_0^2*BETA0pp;
    G.rr(ii)=Amp*alpha_vec(1);
    G.pp(ii)=Amp*alphapp_vec(1);
end


%% Grr
% CxX3=CxX(1:3)
% BETA1_3=BETA1(1:3,1:3);
% BETA0_3=BETA0(1:3);
% inv(diag(CxX)-J_0^2*BETA1)*J_0^2*BETA0
% inv(diag(CxX3)-J_0^2*BETA1_3)*J_0^2*BETA0_3
% inv(diag(CxX(1))-J_0^2*BETA1(1))*J_0^2*BETA0(1)

beta_openloop=@(B0,B1) inv(diag(CxX))*(J_0^2*B1-diag(CxX));
