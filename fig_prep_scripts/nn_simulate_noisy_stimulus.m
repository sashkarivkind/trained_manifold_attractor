function result = nn_simulate_noisy_stimulus(sim,net,xinit,inputmag,epsilonmag,omega,tswitch)
dt=0.1;
tmax=400;

if nargin<4
    inputmag=1e-1;
end
if nargin <5
    epsilonmag=0.5e-1;
end
if nargin<6
    omega=0; %by default - white noise, if frequency set - sine wave at the specified frequency
end
if nargin<7
    tswitch=[80];
end

psi=0;

nsteps=(round(tmax/dt));
u_in=zeros(2,nsteps);
ii=0;
ttt =tswitch;
ii=ii+1;
if omega==0
    u_in(:,round(ttt/dt):end)=[cos(psi);sin(psi)]*...
        (inputmag + epsilonmag*randn(1,size(u_in(:,round(ttt/dt):end),2)));
else
    u_in(:,round(ttt/dt):end)=[cos(psi);sin(psi)]*...
        (inputmag + epsilonmag*sin(omega*dt*(0:(size(u_in(:,round(ttt/dt):end),2)-1))));
    
end
result.noisy=struct;
[result.noisy.xdelta,...
    result.noisy.zdelta,...
    result.noisy.xrecdelta]=fast_conv_to_fp_extended(net,[],...
    struct('xinit',xinit,...
    'ol',0,...
    'tmax',tmax,...
    'dt',dt,...
    'u_in',u_in,...
    'save_neurons',[1:50:50]));


[result.rotation.xdelta_ol_ref,...
    ]=fast_conv_to_fp(...
    net,sim.f_ol+0*u_in(:,end)*ones(1,size(sim.f_ol,2)),...
    struct('ol',1));
result.rotation.zdelta_ol_ref=net.wout'*net.phi(result.rotation.xdelta_ol_ref);

