% cx- correlation of x
% Cr - correaltion of rates
% sigx- sd of x

Amp=1.2;
g_vec=0.1:0.1:2;
M=12; % 6 or 3?....
sim_resolution=180;
ind=0;

% 1. Calculate MF Sampling
% G1p=[];
% G1p1Aprox1=[];
% G1p1Aprox2=[];
for g=g_vec 
    ind=ind+1


psi_vec=(0:(sim_resolution-1))./sim_resolution*2*pi;
net.phi=@(x)erf(x./sqrt(2));
% Solve mean field to get correaltion function
[mu, sigxa,cxa,ca]=...
    solveFixedPointRing_Erf(net.phi,0,g,Amp,1,1,100,psi_vec);
sigx=sqrt(sigxa.^2-Amp^2/2); cx=cxa-Amp.^2.*ca; Cr=cx./g^2;

% Do FFT to get the Fourier
% [cSortTEST,CsimTEST]=compareCorFuncSimMFT(g)

cfu=fft(Cr);
figure(21)
semilogy(real(cfu),'-o')
hold all

% cSort=sort(abs(cfu/length(cfu)),'descend'); % NATOICE FACTOR 2 FROM SASHA SCRIPT WITH COSINE COEFF!
% cSort=sort(2*real(cfu/length(cfu)),'descend'); % NATOICE FACTOR 2 FROM SASHA SCRIPT WITH COSINE COEFF!
% cSort=2*real(cfu(2:2:end)/length(cfu)); 
cSort=2*real(cfu(2:floor(end/2))/length(cfu)); 

% % cSort=CsimTEST';
% % cSort=cSort(1:2:end); % Take off multiplicity due to the fact that C is real and symmetric
% cSort=cSort(2:2:end); % Sasha: WHY did you start from  2 in your MFT???

% semilogy(cSort,'--k');
% legend('sim','MFT')

% Sampling MFT equation
cSortp1=cSort(M+1:M:end);cSortm1=cSort(M-1:M:end-2);
k=1:length(cSortp1);
num=1- (1/cSort(1)).*sum((M.*k-1).*cSortm1-(M.*k+1).*cSortp1);
denum=1+ (1/cSort(1)).*sum((cSortm1-1)+(cSortp1+1));
G1p(2,ind)=num/denum; % I think that here there are issues with numerics- the presicion of the integral for the corr function

num=1- (1/cSort(1)).*((M-1).*cSortm1(1)-(M+1).*cSortp1(1));
denum=1+ (1/cSort(1)).*(cSortm1(1)+cSortp1(1));
G1p1Aprox1(2,ind)=num/denum; % this should give a good approximation of the result
% % G1p1Aprox2(ind)=1-M*cSort(M);

end
figure;%(128)
% plot(g_vec,G1p)
hold all;
plot(g_vec,G1p1Aprox1,'--r')
plot(g_vec,G1p,'--g')
% plot(g_vec,G1p1Aprox2)

% 2. Calculate MF directly from the raw equation

% 3. Calculate G1p by linearizing and look at polar plots

% 4. Calculate G1p by looking at derivative of reconstruction error