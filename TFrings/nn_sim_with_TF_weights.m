clear;
hp=struct;
sim=struct;
net.N=1000;
[hp, net, sim] = prep_network_param(hp, net, sim);
W = readtable('ring_weights_fig103.csv','ReadVariableNames',0);
% rarch = readtable('activity103longDense.csv','ReadVariableNames',0);
W=W{:,:}';
% rarch = rarch{:,:};
[u,s,v]=svd(W);
[~,ii] = sort(diag(s),'descend');

%%
stag=s;
stag(1,1)=0;
stag(2,2)=0;

net.wfb=u(:,1:2)*sqrt(net.N);
net.wout=v(:,1:2)*s(1:2,1:2)/sqrt(net.N);
net.W=u*stag*v';
%%
Q=[3,4];
wfbQ=u(:,Q)*sqrt(net.N);
woutQ=v(:,Q)*s(Q,Q)/sqrt(net.N);
%%
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

%learning output weights
sim.r = net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
regfac = hp.M*hp.alpha_reg*eye(length(pts));
% net.wout = (sim.r(:,pts)/...
%     (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
%%
% net.wout=0*net.wfb;
%obtaining open loop output
sim.z_ol = net.wout'*sim.r;

%simulating closed loop
x_test_rand=fast_conv_to_fp(net,[],struct('xinit',0.05*randn(size(x))));
r_test_rand=net.phi(x_test_rand);
sim.z_test_rand= net.wout'*r_test_rand;

x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
r_test_noi=net.phi(x_test_noi);
sim.z_test_noi= net.wout'*r_test_noi;

%plotting results
figure;
plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
hold on;
plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
%%
[vv,ee]=eig(net.W);
[~,ii]=sort(real(diag(ee)),'descend');

vvr=real(vv(:,ii));

figure;
plot(vvr(:,1)'*x_test_rand,vvr(:,2)'*x_test_rand,'x');
%%
figure;
plot(eig(net.W*diag(net.phip(x(:,23)))),'x')
%%
%%
figure;
plot(eig(net.W*diag(net.phip(x_test_rand(:,10)))),'x')
%%

%%
eee=[];
for uu=1:size(x_test_rand,2)
eee(end+1)=max(real(eig(W*diag(net.phip(x_test_rand(:,uu))))));
end

%%
% eee0=[];
% for uu=1:size(x,2)
% eee0(end+1)=max(real(eig(net.W*diag(net.phip(x(:,uu))))));
% end
% %%


% figure; 
% plot(diag(s),'x')


%% plotting results
figure;
plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
hold on;
zQ=woutQ'*r_test_rand;
plot(zQ(1,:),zQ(2,:),'x', 'linewidth',1);

%% plotting results
figure;
plot(sim.z_ol(1,:),sim.z_ol(2,:),'x', 'linewidth',1);