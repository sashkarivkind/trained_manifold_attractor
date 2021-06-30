clear;
% M_vec=3:2:13;
% A_vec=[1.2]%,1,3,10];
M_vec=60;%[3,5,7,9,11,51];
J_vec=1.5;%g_vec;%[0.2,1.0];%g_vec;%[1.0:0.1:2.1];
phi=@(x)erf(x./sqrt(2));
phip =  @(x) exp(-x.^2/2)/sqrt(pi/2);
% h_vec=[2:2:10];
h_vec=[2:2:4];
omega_vec=[0,logspace(-2,2,40)];
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
        for hh=1:length(h_vec)
%         Amp = A_vec(jj);
        J_0 = J_vec(jj);
        if hh==1
        sfp=0;
        else
        sfp=sfpr{jj,mm,1};    
        end
        [G,debuu{jj,mm,hh},sfpr{jj,mm,hh}]=nn_MFT_Fourier(phi,J_ext,J_0,Amp,K,psi_vec,phip,omega_vec,0,h_vec(hh));
        Grr_rec{jj,mm,hh}=G.rr;
        Gpp_rec{jj,mm,hh}=G.pp;
        Grr0(jj,mm,hh)=Grr_rec{jj,mm,hh}(1);
%         Grr0_unnorm(jj,mm)=debu2(1);
%         lambda_est(jj,mm)=(1-Grr0(jj,mm))*(1- debuu{jj,mm}{2}.kappa);
    end
end
end
% figure; plot(M_vec,Grr0);
%%
figure;
for hh=1:length(h_vec);
subplot(1,3,1); plot(Grr_rec{1,1,hh}); hold on
title('G_r_r')

subplot(1,3,2); plot(abs(Grr_rec{1,1,hh})); hold on
subplot(1,3,3); plot(angle(Grr_rec{1,1,hh})); hold on
end
%%
figure;
for hh=1:length(h_vec);    
subplot(1,3,1); plot(Gpp_rec{1,1,hh}); hold on; 
title('G_{\psi\psi}')
subplot(1,3,2); plot(abs(Gpp_rec{1,1,hh})); hold on
subplot(1,3,3); plot(angle(Gpp_rec{1,1,hh})); hold on
end
% 
%%
figure;
for hh=1:length(h_vec);    
subplot(1,3,1); plot(Gpp_rec{1,1,hh}); hold on; 
title('G_{\psi\psi}')
subplot(1,3,2); loglog(omega_vec,abs(Gpp_rec{1,1,hh})); hold on
subplot(1,3,3); semilogx(omega_vec,angle(Gpp_rec{1,1,hh})); hold on
end

%%
rec_rec=[];
for hh=length(h_vec);
b0=debuu{1,1,hh}.betapp{1}(1,1);
b1=debuu{1,1,hh}.betapp{2}(1,1);
CxX=debuu{1,1,hh}.CxX;
G0=Amp*J_0^2*b0/CxX(1)/(1-J_0^2*b1/CxX(1));
num=Amp*J_0^2*b0/CxX(1);
denum=(1-J_0^2*b1/CxX(1));
rec_rec(end+1,:)=[G0,num,denum];
end
%%
figure;
plot(rec_rec);
%%

phase_games1=1./(debuu{1,1,5}.pp_p1.v1(1,1)+1i*omega_vec);
phase_games2=1./(debuu{1,1,5}.pp_p1.v2(1,1)+1i*omega_vec);
figure; 
plot(angle(phase_games1));
hold on;
plot(angle(phase_games2));
%%

phase_games1=1./(debuu{1,1,1}.pp_p1.v1(1,1)+1i*omega_vec);
phase_games2=1./(debuu{1,1,1}.pp_p1.v2(1,1)+1i*omega_vec);
figure; 
plot(angle(phase_games1));
hold on;
plot(angle(phase_games2));

%
% figure;
% Gpp_ref = Gpp_rec{1,1}(1)./(1+1i*omega_vec);
% subplot(1,3,1); plot(Gpp_rec{1,1},'r');
% hold on; plot(Gpp_ref,'k-');
% 
% subplot(1,3,2); plot(abs(Gpp_rec{1,1}),'r');
% hold on; plot(abs(Gpp_ref),'k-');
% xlabel('\omega'); ylabel('abs(Gpp)');
% subplot(1,3,3); plot(angle(Gpp_rec{1,1}),'r');
% hold on; plot(angle(Gpp_ref),'k-');
% xlabel('\omega'); ylabel('phase(Gpp)');

% legend(num2str(J_vec'))
% xlabel('M');
% ylabel('Grr (\omega=0)');
% figure; plot(M_vec,Grr0_unnorm);
% figure; plot(M_vec,lambdNorm./M_vec);