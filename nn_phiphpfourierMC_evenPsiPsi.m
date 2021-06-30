function [BETA0,BETA1]=phiphpfourierMC_evenPsiPsi(phi,phip,c,A,nmc)
%computing the expected value of integrals
%\int dpsi/2pi \int theta \int c_1..m \int  b_1..m
%phi(theta,psi)phip(theta,0)a_i [or b_i]
%
%Inputs:
% phi - activation function
% phip - activation derivative
% c - relevant coefficients of Fourier
% nmc - number of monte-carlo samples
% A - stimulus amplitude.
% in future planning to support more than one input coefficient
if nargin < 5
    nmc=1e7;
end

tp=rand(nmc,2)*2*pi; % distribution of angles theta and psi
reCo=real(sqrt(real(c(1:2:end)))); %real odd fourier coefficients
ai=randn(nmc,length(reCo))*diag(reCo); %cos coefficient
bi=randn(nmc,length(reCo))*diag(reCo); %sin coefficient
harmonics=tp(:,2)*[1:2:length(c)];
x1_approx=sum((ai.*cos(harmonics))...
    +(bi.*sin(harmonics)),2);

BETA0=zeros(length(reCo),1);
BETA1=zeros(length(reCo));
for ii=1:length(reCo)
    BETA0(ii)=2*mean(...
        sin((2*ii-1)*tp(:,2)).*...
        phi(A*cos(tp(:,1)-tp(:,2))+x1_approx)...%phi integrated over psi
        .*phip(A*cos(tp(:,1))+sum(ai,2))...%at origin we have a sum of all the cosine coeffitients
        .*sin(tp(:,1))...
        );
    for jj=1:length(reCo)
        BETA1(ii,jj)=2*mean(...
            sin((2*ii-1)*tp(:,2)).*... %todo - hard coded only odd harmonics
            phi(A*cos(tp(:,1)-tp(:,2))+x1_approx)...%phi integrated over psi
            .*phip(A*cos(tp(:,1))+sum(ai,2))...%at origin we have a sum of all the cosine coeffitients
            .*bi(:,jj)...
            );
    end
end