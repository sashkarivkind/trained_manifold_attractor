clear;%close all;
%direction - 1 - transversal, 2- parallel
A=1.2;
hp.M=40;


net.g=1.0;
net.N=2000;

net.phi=@(x)erf(x./sqrt(2));
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
    
h_vec=[-20:2:-2,-1:0.1:0.9];
K_vec=1:1:20;
hp.alpha_reg=1e-8;
[hp, net, sim] = prep_network_param(hp,net,struct);

for hh=1:length(h_vec)
    h=h_vec(hh);

    
    sim.f_ol=A*[cos(sim.psi); (1-h)*sin(sim.psi)];
    
    x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
    sim.x=x;
    %learning output weights
    sim.r = net.phi(x);
    pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
    regfac = net.N*hp.M*hp.alpha_reg*eye(length(pts));
    net.wout = (sim.r(:,pts)/...
        (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
    wo_rec{hh}=net.wout;
    %obtaining open loop output
    sim.z_ol = net.wout'*sim.r;
    
    %simulating closed loop
    % x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
    % r_test_rand=net.phi(x_test_rand);
    % sim.z_test_rand= net.wout'*r_test_rand;
    %
    % x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
    % r_test_noi=net.phi(x_test_noi);
    % sim.z_test_noi= net.wout'*r_test_noi;
    %
    % %plotting results
    % figure;
    % plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
    % hold on;
    % plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
    % plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
    % plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
    % plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
    %
    %
    pts=1:hp.sim_resolution;
    pt=1;
    % ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
    % ee_cl1=eig((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt))));
    % ee_cl2=eig((net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(x(:,pt))));
    % % ee_clg0=eig((0*net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
    %
    % figure;
    % plot(ee_cl-1,'*');
    % hold on;
    %
    % plot(real(ee_cl1(1:5)-1),imag(ee_cl1(1:5)-1),'o');
    % plot(real(ee_cl2(1:2)-1),imag(ee_cl2(1:2)-1),'s');
    
    legend( 'full','transversal','tangent');
    %%
    K_vec=1:1:20;
    for kk=1:length(K_vec)
        K=K_vec(kk);
        o{hh,kk} = semi_empirical_spectrum(net,sim,K,pts,pt);
    end
hh    
end
    % hold on;plot(y,'x','MarkerSize',15)
    % figure;imagesc(MM)
    % figure;imagesc(inv(lhs)*B1);title('beta1')
    
