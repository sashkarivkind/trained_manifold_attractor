function [result]=nn_plot_compare_to_semiemp(a,h,A,g,N,M,seed,UniformFlag,PlotFlag,PlotEigFlag,K_vec)
% 22/01/21
% Function runs the simulaiton of te learned network and is used later to
% plot te tuning curves
% UniformFlag - for uniform distribution of points on te manifold
% PlotFlag    - to be sure that manifold selceted is stable
% seed        - to be able to replicate

if nargin <11
    K_vec=[2:2:30];
end
rng(seed)
result=struct;
tmax=1000;
dt=1; b=1;rot=0;
net.g=0; net.N=N;
net.phi  = @(x)erf(x./sqrt(2)); net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
net.phip = @(x) exp(-x.^2/2)/sqrt(pi/2);

hp=struct;
sim=struct;
hp.A=A;

[hp, net, sim] = prep_network_param(hp, net,sim);
hp.M=M;
if UniformFlag
    tmp_f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
    sim.f_ol=curvspace(tmp_f_ol',hp.sim_resolution);
    sim.f_ol=sim.f_ol';
else
    sim.f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
end

x0 = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
trained_net = update_net_g(net,g);
x = fast_conv_to_fp(trained_net,sim.f_ol,struct('ol',1));
sim.r = trained_net.phi(x);
pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
trained_net.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
result.net=trained_net;
if PlotEigFlag
    % result.ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
    result.ee_cl1=eig((trained_net.W+trained_net.wfb(:,1)*trained_net.wout(:,1)')*diag(trained_net.phip(x(:,1)))-eye(trained_net.N));
    result.ee_cl2=eig((trained_net.W+trained_net.wfb(:,2)*trained_net.wout(:,2)')*diag(trained_net.phip(x(:,1)))-eye(trained_net.N));
    figure
    plot(result.ee_cl1,'.g')
    hold on
    plot(result.ee_cl2,'.r')
    %obtaining the semi-MF eigenvalues
    ptsMF=1:1:length(sim.f_ol);
    
    for kk=1:length(K_vec)
        for dd=1:2
        
            K=K_vec(kk);
            sim.x=x;
            semi=semi_empirical_spectrum(trained_net,sim,K,1:160,1,dd);
            
            result.y{kk,dd}=eig(semi.MM)-1;
        end
        semi_all=semi_empirical_spectrum(trained_net,sim,K,1:160,1);
            result.y_all{kk}=eig(semi_all.MM)-1;

    end
    sig_0=1;[sigma]=solveFixedPointMuSigRing_Erf(trained_net.phi,trained_net.g,hp.A,sig_0);
    result.bulkMF=bulk_fun(trained_net.g,sigma,hp.A,trained_net.phip);
    %         rplot=1+result.bulkMF;xplot=-1;yplot=0; th = 0:pi/60:2*pi;xunit = rplot* cos(th) + xplot; yunit = rplot * sin(th) + yplot;
    %         plot(xunit, yunit,'--k');
end

result.with_learning.z_ol = trained_net.wout'*sim.r;
result.sim = sim;
results.K_vec=K_vec;
if PlotFlag
    figure
    %     x_test_rand=fast_conv_to_fp(trained_net,[],struct('xinit',5*randn(size(x))));
    %     r_test_rand=trained_net.phi(x_test_rand);
    %     sim.z_test_rand= trained_net.wout'*r_test_rand;
    
    x_test_noi=fast_conv_to_fp(trained_net,[],struct('xinit',x+1e-3*randn(size(x))));
    r_test_noi=trained_net.phi(x_test_noi);
    sim.z_test_noi= trained_net.wout'*r_test_noi;
    
    %plotting results
    plot(result.with_learning.z_ol(1,:),result.with_learning.z_ol(2,:),'-r');
    hold on;
    %     plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
    %     plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
    
    %     plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'og', 'linewidth',1);
    colorVec = hsv(size(sim.z_test_noi,2));
    for ll=1:size(sim.z_test_noi,2)
        plot(sim.z_test_noi(1,ll),sim.z_test_noi(2,ll),'o', 'Color',colorVec(ll,:));
    end
    %     box off; axis square;
    %         plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'xb', 'linewidth',1);
    
    xlim([-2 2]); ylim([-2 2])
    title(['A= ' num2str(A),'g= ' num2str(g),'a= ' num2str(a),'h= ' num2str(h) ])
    
end


