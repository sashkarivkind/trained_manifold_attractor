run_new_net = 1;
if run_new_net
    clear;
end
    
    
    %setting hyper parameters
    K=10;
    pt=1;
    hp.alpha_reg = 0.0;
    hp.b=1;
    hp.M=20;
    hp.sim_resolution=160; %points to simulate on the ring
    hp.omega_vec = logspace(-2,4,40);
    rng(1);
    A_vec=[0.4:0.1:1.5];
    g_vec=[0.01, 0.1:0.2:1.3];

    %setting network parameters
    for aa=1:length(A_vec)
    for gg=1:length(g_vec)
        net.g=g_vec(gg);
        hp.A=A_vec(aa);

        net.N=1000;
%         net.win=randn(net.N,1);
        % net.phi=@(x) tanh(x); %added 1 to break symmetry
        % net.phip =@(x) 1-tanh(x).^2;
        
        net.phi=@(x)erf(x./sqrt(2));
        net.phip = @(x) exp(-x.^2/2)/sqrt(pi/2);
        
        net.W=net.g*randn(net.N)/sqrt(net.N);
        theta=2*pi*(0:net.N-1)/net.N;
        net.wfb=randn(net.N,2);%[cos(theta'),sin(theta')];
        
        sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate
        sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        
        %learning output weights
        r = net.phi(x);
        
        pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
        regfac = hp.M*hp.alpha_reg*eye(length(pts));
        net.wout = (r(:,pts)/(r(:,pts)'*r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
        
        %obtaining open loop output
        sim.x=x;
        sim.r=r;
        sim.z_ol = net.wout'*r;
        MlinR=(net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,1)));
%        ee=sort(MlinR,'descend') 
        lambdaR(aa,gg)=max(real(eig(MlinR)));
        sigmax(aa,gg)=sqrt(mean(x(:,1).^2));
        pts_for_semi=1:length(sim.psi);
        semi=semi_empirical_spectrum(net,sim,K,pts_for_semi,pt,1);
        semi2=semi_empirical_spectrum(net,sim,2,pts_for_semi,pt,1);
        lambdaRsemi(aa,gg)=max(real(eig(semi.MM)));
        lambdaRsemi2(aa,gg)=max(real(eig(semi2.MM)));

    disp([aa,gg]);
    end
    end