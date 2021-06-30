% cx- correlation of x
% Cr - correaltion of rates
% sigx- sd of x

Amp=1.2;
g_vec=[1.2:0.2:2];
M=6; % 6 or 3?....
sim_resolution=160;
ind=0;
Kmax=10;
%Calculate MF Sampling
G1p=[];
G1p1Aprox1=[];
G1p1Aprox2=[];
for g=g_vec 
    ind=ind+1


psi_vec=(0:(sim_resolution-1))./sim_resolution*2*pi;
net.phi=@(x)erf(x./sqrt(2));
% Solve mean field to get correaltion function
[mu, sigxa,cxa,ca]=...
    solveFixedPointRing_Erf(net.phi,0,g,Amp,1,1,100,psi_vec);
sigx=sqrt(sigxa.^2-Amp^2/2); cx=cxa-Amp.^2.*ca; Cr=cx./g^2;

cfu=fft(Cr);
cfudif=fft(diff(Cr)./(psi_vec(2)-psi_vec(1)));
G1p1Aprox0(ind)=abs(-(cfudif(2))/(cfu(2)))

figure(21)
semilogy(real(cfu),'-o')
hold all

cSort=2*real(cfu(2:2:end)/length(cfu)); 
% Sampling MFT equation
cSortp1=cSort(M+1:M:end);cSortm1=cSort(M-1:M:end-2);
k=1:length(cSortp1);
num=1- (1/cSort(1)).*((M-1).*cSortm1(1)-(M+1).*cSortp1(1));
denum=1+ (1/cSort(1)).*(cSortm1(1)+cSortp1(1));
G1p1Aprox1(ind)=num/denum 
% G1p1Aprox2(ind)=1-M*cSort(M)
% this should give a good approximation of the result

k=1:length(cSortp1);
num=1- (1/cSort(1)).*sum((M.*k(1:Kmax)-1).*cSortm1(1:Kmax)-(M.*k(1:Kmax)+1).*cSortp1(1:Kmax));
denum=1+ (1/cSort(1)).*sum((cSortm1(1:Kmax)+(cSortp1(1:Kmax))));
G1p(ind)=num/denum 
% I think that here there are issues with numerics- the presicion of the integral for the corr function


end
figure(3)

hold all;
plot(g_vec,G1p1Aprox0,'--b')
plot(g_vec,G1p1Aprox1,'--g')
plot(g_vec,G1p,'-g')
