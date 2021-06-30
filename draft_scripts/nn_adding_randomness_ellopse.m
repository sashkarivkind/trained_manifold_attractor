clear;
result=struct;
tmax=1000;
dt=1;
h=0.1;
A=1.2;
net.g=1.0;
net.N=1000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);

[hp, net, sim] = prep_network_param(struct, net,struct);
hp.alpha_reg=1e-8;
sim.f_ol=A*[cos(sim.psi); (1-h)*sin(sim.psi)];


%% learning ring with g
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

sim.r = net.phi(x);
result.with_learning={};
M_vec=[5,20];
for MM=1:length(M_vec)
hp.M=M_vec(MM);

pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = net.N*hp.M*hp.alpha_reg*eye(length(pts));
net.wout = (sim.r(:,pts)/...
    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution

result.with_learning{MM}.z_ol = net.wout'*sim.r;
[x_cl,result.with_learning{MM}.z_cl]=fast_conv_to_fp_extended(net,[],...
struct('xinit',x,...
'ol',0,...
'tmax',tmax,...
'dt',dt));
end

figure;
plot(sim.psi,phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi);
hold on;
plot(sim.psi,phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi);
% legend('classic','small niose','large noise','M=5','M=20');

figure;plot(angle(result.with_learning{1}.z_cl)');
figure;plot(angle(result.with_learning{2}.z_cl)');
figure; plot(result.with_learning{2}.z_cl(:,end),'x'),
