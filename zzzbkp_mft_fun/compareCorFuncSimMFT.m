function [cSort,Csim]=compareCorFuncSimMFT(g)
sim_resolution=160; %todo: change in sim struct
Amp=1.2; %todo: change in sim struct
% Issues: 1. the MF solution gives the spectrum upt to ~ 1e-7
%         2. the simulations deviated from MF for large g (e.g. differences in pairs of e.v.)
%            -probably due to finite N corrections (sim_res is not large enough). 

% 1. Run simulation
net.phi=@(x)erf(x./sqrt(2));
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
net.g=g;
net.N=1000;
% hp.sim_resolution=sim_resolution;
[hp, net, sim] = prep_network_param(struct,net,struct);
net.N
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
%learning output weights
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
%obtaining open loop output
sim.z_ol = net.wout'*sim.r;
result.C=(sim.r'*sim.r)/net.N/hp.sim_resolution;
%plotting results
figure;
Csim=flip(eig(result.C));
semilogy(Csim,'b');
hold on;

% 2. Calculate MF
psi_vec=(0:(sim_resolution-1))./sim_resolution*2*pi;
net.phi=@(x)erf(x./sqrt(2));
% Solve mean field to get correaltion function
[mu, sigxa,cxa,ca]=...
    solveFixedPointRing_Erf(net.phi,0,g,Amp,1,1,100,psi_vec);
sigx=sqrt(sigxa.^2-Amp^2/2); cx=cxa-Amp.^2.*ca; Cr=cx./g^2;
cfu=fft(Cr);
cSort=sort(abs(cfu/length(cfu)),'descend'); % NATOICE FACTOR 2 FROM SASHA SCRIPT WITH COSINE COEFF!
semilogy(cSort,'--k');
legend('sim','MFT')

end