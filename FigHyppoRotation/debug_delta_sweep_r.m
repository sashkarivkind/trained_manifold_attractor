r_vec=0.9:0.1:1.5;
clear z_epsi0 z_epsi
for rr=1:length(r_vec)
    
f_ol_tag=r_vec(rr)*[cos(sim.psi);sin(sim.psi)];

x_epsi0 = fast_conv_to_fp(net,...
    f_ol_tag,...
    struct('ol_with_fixed_input',1,...
    'u',0*u_in(:,1000)*ones(1,size(sim.f_ol,2))));
z_epsi0{rr} = net.wout'*net.phi(x_epsi0);

x_epsi = fast_conv_to_fp(net,...
    f_ol_tag,...
    struct('ol_with_fixed_input',1,...
    'u',u_in(:,1000)*ones(1,size(sim.f_ol,2))));
z_epsi{rr} = net.wout'*net.phi(x_epsi);

end

%%
figure;

for rr=1:length(r_vec)
plot(180/pi*sim.psi,angle(z_epsi0{rr}(1,:)+1i*z_epsi0{rr}(2,:))-sim.psi,'x-')
hold on
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1);
plot(180/pi*sim.psi,angle(z_epsi{rr}(1,:)+1i*z_epsi{rr}(2,:))-sim.psi,'o-')
rr
end
grid