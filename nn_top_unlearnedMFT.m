% M_vec=3:2:13;
% A_vec=[1.2]%,1,3,10];
M_vec=160;%[3,5,7,9,11,51];
J_vec=g_vec;%[0.2,1.0];%g_vec;%[1.0:0.1:2.1];
phi=@(x)erf(x./sqrt(2));
phip =  @(x) exp(-x.^2/2)/sqrt(pi/2);
omega_vec=[0,1];%logspace(-2,2,40)];
Amp = 1.2;
J_ext=0;
J_0=0;
K=100;
Gtt_rec={};
debuu={};
sfpr={};
Grr0=[];
Grr0_unnorm=[];
lambdNorm=[];
lambda_est = [];
for mm=1:length(M_vec)
    Mtrain = M_vec(mm);
    psi_vec=2*pi*(0:Mtrain-1)/Mtrain;
    for jj=1:length(J_vec)
%         Amp = A_vec(jj);
        J_0 = J_vec(jj);
        [G,debuu{jj,mm},sfpr{jj,mm}]=nn_MFT_Fourier(phi,J_ext,J_0,Amp,K,psi_vec,phip,omega_vec,0,10);
        Grr_rec{jj,mm}=G.rr;
        Gpp_rec{jj,mm}=G.pp;
        Grr0(jj,mm)=Grr_rec{jj,mm}(1);
%         Grr0_unnorm(jj,mm)=debu2(1);
%         lambda_est(jj,mm)=(1-Grr0(jj,mm))*(1- debuu{jj,mm}{2}.kappa);
    end
end

% figure; plot(M_vec,Grr0);
% figure;
% subplot(1,3,1); plot(Grr_rec{1,1},'r');
% subplot(1,3,2); plot(abs(Grr_rec{1,1}),'r');
% subplot(1,3,3); plot(angle(Grr_rec{1,1}),'r');

figure;
Gpp_ref = Gpp_rec{1,1}(1)./(1+1i*omega_vec);
subplot(1,3,1); plot(Gpp_rec{1,1},'r');
hold on; plot(Gpp_ref,'k-');

subplot(1,3,2); plot(abs(Gpp_rec{1,1}),'r');
hold on; plot(abs(Gpp_ref),'k-');
xlabel('\omega'); ylabel('abs(Gpp)');
subplot(1,3,3); plot(angle(Gpp_rec{1,1}),'r');
hold on; plot(angle(Gpp_ref),'k-');
xlabel('\omega'); ylabel('phase(Gpp)');

% legend(num2str(J_vec'))
% xlabel('M');
% ylabel('Grr (\omega=0)');
% figure; plot(M_vec,Grr0_unnorm);
% figure; plot(M_vec,lambdNorm./M_vec);