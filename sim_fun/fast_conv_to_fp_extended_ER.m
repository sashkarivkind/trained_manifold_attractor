function [xout,z,xrec]=fast_conv_to_fp_extended_ER(net,f_ol,opt)

nom.tmax=1000;
nom.xinit = nan;
nom.save_neurons=[];
nom.ol=0;
nom.dt=1;
nom.u =zeros(size(net.W,1),1); %legacy, to be removed after ensuring we don't use it anywhere
nom.u_in = 0;
opt=nom_opt_assigner(opt,nom);
clear nom;

%attaching zero input weights to a network if it does not have input wieghts
nom_net.win=0;
net = nom_opt_assigner(net,nom_net);

nsteps=round(opt.tmax/opt.dt);
opt.u_in=opt.u_in.*ones(size(opt.u_in,1),nsteps);
if isnan(opt.xinit)
    xout=randn(net.N,nsteps);
else
    xout = opt.xinit;
end

xrec=zeros(length(opt.save_neurons),nsteps);
z=zeros(size(xout,2),nsteps);
for t=1:nsteps
    z(:,t)=[net.wout(:,1)'*net.phi(xout) + 1i*net.wout(:,2)'*net.phi(xout)].';
    if opt.ol
        xout=xout+opt.dt*(...
            -xout+net.W*net.phi(xout)...
            +net.wfb*f_ol...
            +net.win*opt.u_in(t)...
            +opt.u*ones(1,size(f_ol,2))+net.BalanceInp);
    else
        xout=xout+opt.dt*(...
            -xout...
            +(net.W+net.wfb*net.wout')*net.phi(xout)...
            +net.win*opt.u_in(:,t)...
            +opt.u*ones(1,size(xout,2))+net.BalanceInp);
    end
    xrec(:,t)=xout(opt.save_neurons);
end

