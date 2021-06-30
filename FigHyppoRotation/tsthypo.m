clear;
for seed=1;%8:10;
    seed
run_new_net = 1;
zdelta_rec={};
% g_vec=0:0.2:2;
g_vec=[0.2];%0:0.1:2;
input_mode = 'hyppo';
t_conv_vec=[];
for gg=1:length(g_vec)
   tic;
    gg
    if 1 %
        rng(seed);
        hp=struct;
        net=struct;
        hp.M=40;
        A=0.9;
        a=4;b=1;
        h=0.7;
        rot=0;
        sim=struct;

        net.g=g_vec(gg);
        net.N=1000;
        
        [hp, net, sim] = prep_network_param(hp, net, sim);
        
        sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));

        
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        %learning output weights
        sim.r = net.phi(x);
        pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
        regfac = hp.M*hp.alpha_reg*eye(length(pts));
        net.wout = (sim.r(:,pts)/...
            (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
        
        %obtaining open loop output
        sim.z_ol = net.wout'*sim.r;
        
%         %simulating closed loop
%         x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
%         r_test_rand=net.phi(x_test_rand);
%         sim.z_test_rand= net.wout'*r_test_rand;
%         
%         x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
%         r_test_noi=net.phi(x_test_noi);
%         sim.z_test_noi= net.wout'*r_test_noi;
        
        %plotting results
%         figure;
%         plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
%         hold on;
%         plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
%         plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
%         plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
%         plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
     %%   
    end
    net.win=net.wfb;
    
    dt=1;
    tmax=6000;
%     epsilonmag=0.01;
    epsilonmag=0.001;
    tswitch=10*[50, 250, 450];
    pulse_width = 100./dt*[10,10,10];
    psi_vec=[75,130, 30]*pi/180;
%     tswitch=10*[50, 250];
%     pulse_width = 1000*[100,100];
%     psi_vec=[135 120]*pi/180;
    nsteps=(round(tmax/dt));
    %spikes
    u_in=zeros(2,nsteps);
    ii=0;
    for ttt =tswitch
        ii=ii+1;
        psi=psi_vec(ii);
        uu1=round(ttt/dt);
        uu2=min(size(u_in,2),uu1+pulse_width(ii));
        %%
        if strcmp(input_mode,'sincos')
        u_in(:,uu1:uu2)=epsilonmag*[cos(psi);sin(psi)]*...
            ones(1,size(u_in(:,uu1:uu2),2));
        elseif  strcmp(input_mode,'hyppo')
            u_in(:,uu1:uu2)=...
                epsilonmag*1/(a-b)*...
                [(a-b)*cos(psi-rot)+h*cos((a-b)/b*(psi-rot));...
                (a-b)*sin(psi-rot)-h*sin((a-b)/b*(psi-rot))]./(1+h/(a-b))...
                *ones(1,size(u_in(:,uu1:uu2),2));
        else
            error
        end
        %%
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
        'save_neurons',[1:50:net.N]));
%%
T=dt:dt:tmax;
psi_hat=psi_decoded([real(zdelta);imag(zdelta)],sim.f_ol,sim.psi);    
figure;
%     subplot(2,1,1);
    plot(T,angle(zdelta));
    hold on;
    plot(T,psi_hat);
    hold on;
    plot(T,abs(zdelta))
    plot(T,180/pi*abs(u_in(1,:)+1i*u_in(2,:)),'p');
%     subplot(2,1,2);
%     plot(T,xrecdelta');
%% 
figure;
%     subplot(2,1,1);
    plot(T,180/pi*psi_hat);
    hold on;
    plot(T,abs(zdelta))
%     subplot(2,1,2);
%     plot(T,xrecdelta');
%%
    q=0.01;
    t_conv=abs(angle(zdelta)-angle(zdelta(end)))<q*abs(angle(zdelta(1))-angle(zdelta(end)));
%     t_conv_vec(end+1)=min(find(t_conv))*dt-tswitch;
    zdelta_rec{gg}=zdelta;
end
net.W=[];
net.wout=[];
net.win=[];
save(['tstResultsN',num2str(net.N),'seed',num2str(seed)]);

toc;
end
%%
figure
plot(sim.f_ol(1,:),sim.f_ol(2,:),'or')
hold on
plot(zdelta,'.k')

% save('working_state_input4k.mat')