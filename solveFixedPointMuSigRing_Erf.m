function [sig_out]=solveFixedPointMuSigRing_Erf(phi,J_0,A,sig_0)
% Ver1: 13/01/21 
% Function solves the fixed point statistics of a ring atracttor 
% with erf: r=phi(x)=erf(x/sqrt(2))=int_0^xdye^(-y^2/2)1/sqrt(2pi)? 
%
%
% Inputs: 
%       phi             Transfer function
%       J_0             Synaptic strength
%       A               Amplitude of the bump
%       sig_0      initial conditions for sigma
% Outputs
%       sig_out         sd of the field x1

D=@(z)exp(-((z.^2)/2))./(sqrt(2*pi));
dtheta=1./(2*pi);
options=optimset('Display','off');% Tae off if ou want to see display!
%Solve mu and sigma. sigma is sigx (of x1)
a0=sig_0;
y=fsolve(@F,a0,options);
sig_out=y(1);


    function y=F(a)
        sig=a(1);
        g2=@(theta,z) phi(sig.*z+A*cos(theta)).^2.*D(z).*dtheta;
        y(1)=J_0^2*integral2(g2,0,2*pi,-Inf,Inf)-(sig.^2);
    end

end





