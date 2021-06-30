function xout=fast_conv_to_fp(net,f_ol,opt)

nom.tmax=1000;
nom.xinit = nan;
nom.u =zeros(size(net.W,1),1); 
nom.ol=0;
nom.ol_with_fixed_input=0;
nom.cl_with_fixed_input=0;
opt=nom_opt_assigner(opt,nom);
clear nom;

if isnan(opt.xinit)
xout=randn(net.N,size(f_ol,2));
else
xout = opt.xinit;
end

for t=1:opt.tmax
    if opt.ol
            xout=net.W*net.phi(xout)+net.wfb*f_ol+opt.u*ones(1,size(f_ol,2));
    elseif opt.ol_with_fixed_input %this case is defined as a separate option to ensure backward compatibility
            xout=net.W*net.phi(xout)+net.wfb*f_ol+net.win*opt.u;
    elseif opt.cl_with_fixed_input %this case is defined as a separate option to ensure backward compatibility
            xout=(net.W+net.wfb*net.wout')*net.phi(xout)+net.win*opt.u;
    else
        xout=(net.W+net.wfb*net.wout')*net.phi(xout)+opt.u*ones(1,size(xout,2));
    end
end