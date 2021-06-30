function result = nn_simulate_tracking(sim,net,xinit,psi,epsilonmag,tswitch)
dt=0.1;
tmax=400;

if nargin<4
psi_vec=[40]*pi/180;
end
if nargin <5
epsilonmag=3e-1;
end
if nargin<6
tswitch=[80];
end
nsteps=(round(tmax/dt));
u_in=zeros(2,nsteps);
ii=0;
ttt =tswitch;
ii=ii+1;
% u_in(:,round(ttt/dt):end)=epsilonmag*[cos(psi);sin(psi)]*...
%                 ones(1,size(u_in(:,round(ttt/dt):end),2));

 u_in=epsilonmag*[cos(psi);sin(psi)];


result.rotation=struct;
[result.rotation.xdelta,...
    result.rotation.zdelta,...
    result.rotation.xrecdelta]=fast_conv_to_fp_extended(net,[],...
    struct('xinit',xinit,...
    'ol',0,...
    'tmax',tmax,...
    'dt',dt,...
    'u_in',u_in,...
    'save_neurons',[1:50:50]));

[result.rotation.xdelta_ol,...
    ]=fast_conv_to_fp(...
    net,sim.f_ol+u_in(:,end)*ones(1,size(sim.f_ol,2)),...
    struct('ol',1));
result.rotation.zdelta_ol=net.wout'*net.phi(result.rotation.xdelta_ol);

[result.rotation.xdelta_ol_ref,...
    ]=fast_conv_to_fp(...
    net,sim.f_ol+0*u_in(:,end)*ones(1,size(sim.f_ol,2)),...
    struct('ol',1));
result.rotation.zdelta_ol_ref=net.wout'*net.phi(result.rotation.xdelta_ol_ref);

