clear
% - Attempt to fix the "knee" shape at exponential decay
% cx- correlation of x
% Cr - correaltion of rates
% sigx- sd of x

Amp=1.2;
% g_vec=[0.01 0.1:0.1:2.2];
g_vec=[0.2,0.5,1.,1.6,1.8];
M_vec=[4:2:14];
% sim_resolution=120;
ind=0;



% 1. Calculate MF Sampling
% G1p=[];
% G1p1Aprox1=[];
% G1p1Aprox2=[];
for g=g_vec 
    ind=ind+1


% psi_vec=(0:(sim_resolution-1))./sim_resolution*2*pi;
% net.phi=@(x)erf(x./sqrt(2));
% Solve mean field to get correaltion function

load(['../mft_fun/data/save_mft_corr_',strrep(num2str(g),'.','p')]...
    ,'mu', 'sigxa','cxa','ca');
% [mu, sigxa,cxa,ca]=...
%     solveFixedPointRing_Erf(net.phi,0,g,Amp,1,1,100,psi_vec);
sigx=sqrt(sigxa.^2-Amp^2/2); cx=cxa-Amp.^2.*ca; Cr=cx./g^2;

% Do FFT to get the Fourier
% [cSortTEST,CsimTEST]=compareCorFuncSimMFT(g)

cfu=fft(Cr);

cSort=2*real(cfu(2:floor(end/2))/length(cfu)); 

for mm=1:length(M_vec)
M=M_vec(mm);% Sampling MFT equation
cSortp1=cSort(M+1:M:end);cSortm1=cSort(M-1:M:end-2);
k=1:length(cSortp1);
num=1- (1/cSort(1)).*sum((M.*k-1).*cSortm1-(M.*k+1).*cSortp1);
denum=1+ (1/cSort(1)).*sum((cSortm1-1)+(cSortp1+1));
G1p(mm,ind)=num/denum; % I think that here there are issues with numerics- the presicion of the integral for the corr function

num=1- (1/cSort(1)).*((M-1).*cSortm1(1)-(M+1).*cSortp1(1));
denum=1+ (1/cSort(1)).*(cSortm1(1)+cSortp1(1));
G1p1Aprox1(mm,ind)=num/denum; % this should give a good approximation of the result

end
%%
end

%recovering open loop pole
ol_pole_vec=zeros(size(g_vec));
beta_info=load('../mft_fun/data/beta_graph_long.mat');
for gg=1:length(g_vec)
    ii=find(abs(beta_info.J_vec-g_vec(gg))<1e-10);
    if length(ii)~=1
        error('pole info not found!');
    end
    ol_pole_vec(gg)=beta_info.debuu{ii}.pp_p1.v1(1);
end

lambdas=(G1p-1)*diag(ol_pole_vec);
drift_expect=diag(1./M_vec)*(1-G1p)*diag(ol_pole_vec);

% figure;%(128)
% % plot(g_vec,G1p)
% hold all;
% plot(g_vec,G1p1Aprox1,'--r')
% plot(g_vec,G1p,'--g')

% figure;
% plot(M_vec,G1p')
% 
% 
% figure;
% semilogy(M_vec,1-G1p')

figure(28);
% semilogy(M_vec,-lambdas');
   set(gca,'ColorOrderIndex',1);
semilogy(M_vec,(abs(-lambdas')));
hold all
% plot(M_vec,lambdas');
figure(5);
% semilogy(M_vec,-lambdas');
   set(gca,'ColorOrderIndex',1);
semilogy(M_vec,(drift_expect),'--');
hold all

% plot(g_vec,G1p1Aprox2)

% 2. Calculate MF directly from the raw equation

% 3. Calculate G1p by linearizing and look at polar plots

% 4. Calculate G1p by looking at derivative of reconstruction error