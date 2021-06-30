run_new_net = 1;
if run_new_net
    clear;
end
    
    
    %setting hyper parameters
    hp.alpha_reg = 0.0;
    hp.b=1;
    hp.M=20;
    hp.A=1.2;
    hp.sim_resolution=160; %points to simulate on the ring
    hp.omega_vec = logspace(-2,4,40);
    rng(1);
    %setting network parameters
    for g=[0.01  1.5]
        net.g=g
        net.N=1000;
%         net.win=randn(net.N,1);
        % net.phi=@(x) tanh(x); %added 1 to break symmetry
        % net.phip =@(x) 1-tanh(x).^2;
        
        net.phi=@(x)erf(x./sqrt(2));
        net.phip = @(x) exp(-x.^2/2)/sqrt(pi/2);
        
        net.W=net.g*randn(net.N)/sqrt(net.N);
        theta=2*pi*(0:net.N-1)/net.N;
        net.wfb=[cos(theta'),sin(theta')];
        
        sim.psi=2*pi*(0:(hp.sim_resolution-1))/hp.sim_resolution; %angles to simulate
        sim.f_ol=hp.A*[cos(sim.psi);hp.b*sin(sim.psi)];
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        
        %learning output weights
        r = net.phi(x);
        
        pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
        regfac = hp.M*hp.alpha_reg*eye(length(pts));
        net.wout = (r(:,pts)/(r(:,pts)'*r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
        
        %obtaining open loop output
        sim.z_ol = net.wout'*r;
        
        %simulating closed loop
        x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
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
        
        
        net.win=net.wfb(:,1);
        % net.win=net.wfb;
        %randn(net.N,1);
        % net.win=randn(net.N,1);
        dt=0.1;
        tmax=50;
        deltamag=.01;
        tspikes=[10:40:tmax  ];
        nsteps=(round(tmax/dt));
        %spikes
        u_in=zeros(1,nsteps);
        for tspike =tspikes
            u_in(round(tspike/dt))=deltamag/dt;
        end
        %tracking
        % dpsi_dt=0.1;
        % t_vec=dt:dt:tmax;
        % u_in=[-sin(t_vec*dpsi_dt) ;cos(t_vec*dpsi_dt) ];
        [xdelta,zdelta,xrecdelta]=fast_conv_to_fp_extended(net,[],...
            struct('xinit',x(:,1),...
            'ol',0,...
            'tmax',tmax,...
            'dt',dt,...
            'u_in',u_in,...
            'save_neurons',[1:50:1000]));
        
        % figure;
        % subplot(2,1,1);
        % plot(angle(zdelta));
        % hold on;
        % plot(abs(zdelta))
        % subplot(2,1,2);
        % plot(xrecdelta');
        net
        figure(221)
        plot(abs(zdelta))
        hold on
    end
    %%
    subplot 211
    xlim ([0 1000])
    ylim ([-0.001 1.5])
    box off
    subplot 212
    xlim ([0 1000])
    ylim ([-3 3])
    box off
