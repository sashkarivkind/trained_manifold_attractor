
UniformFlag=0;
seed=1;
rng(seed);
N=4000;M=40;
g_vec=[1 1.2 0.6];
a_vec=[4 6 2];
h_vec=[0.6 0.9  0.8];
for ii=1:3
    ii
    g=g_vec(ii);h=h_vec(ii);a=a_vec(ii);
    % N=4000;M=40;g=1; A=1.2; h=0.6; a=4;
    % N=4000;M=40;g=1.2; A=1.2; h=0.9; a=6;
%    g=0.6; A=1.2; h=0.8; a=2;
    
    result=struct;
    tmax=1000;
    dt=1; b=1;rot=0;
    trained_net.g=g; trained_net.N=N;
    trained_net.phi  = @(x)erf(x./sqrt(2)); trained_net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
    trained_net.phip = @(x) exp(-x.^2/2)/sqrt(pi/2);
    
    hp=struct;
    sim=struct;
    hp.sim_resolution=160;
    hp.A=A;
    hp.M=M;
    [hp, trained_net, sim] = prep_network_param(hp, trained_net,sim);
    
    if UniformFlag
        tmp_f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
        sim.f_ol=curvspace(tmp_f_ol',hp.sim_resolution);
        sim.f_ol=sim.f_ol';
    else
        sim.f_ol=hp.A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
    end
    
    % trained_net = update_net_g(net,g);
    x = fast_conv_to_fp(trained_net,sim.f_ol,struct('ol',1));
    sim.r = trained_net.phi(x);
    pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
    trained_net.wout=lms_weights(sim.r(:,pts),sim.f_ol(:,pts));
    
    ptind=0;
    for pt=1:1:hp.sim_resolution
        tic;
        pt
        ptind=ptind+1;
        S=(trained_net.W+trained_net.wfb*trained_net.wout')*diag(trained_net.phip(x(:,pt)))-eye(trained_net.N);
        result.ee_cl(ptind,:)=eigs(S,2,'largestreal');
        toc;
    end
    sig_0=1;[sigma]=solveFixedPointMuSigRing_Erf(trained_net.phi,trained_net.g,hp.A,sig_0);
    result.bulkMF=bulk_fun(trained_net.g,sigma,hp.A,trained_net.phip)
    
    %
    jump=1;
    EigNumber=2;
    figure
    plot(1:jump:hp.sim_resolution,result.ee_cl(:,1))
    hold on
    plot(1:jump:hp.sim_resolution,result.ee_cl(:,2))
    plot(1:jump:hp.sim_resolution,result.bulkMF.*ones(length(1:jump:hp.sim_resolution),1))
    
    figure
    scatter(sim.f_ol(1,1:jump:end),sim.f_ol(2,1:jump:end),30,result.ee_cl(:,EigNumber),'o')
    colorbar
    xlim([-1.5 1.5]);ylim([-1.5 1.5])
    set(gcf,'Position',[ 903   775   217   173])
    hAx=gca
    colormap(hsv)
    save(['ExamplesMaxEVAlongmanifold',num2str(trained_net.N),...
                'a',strrep(num2str(a),'.','p'),...
                'h',strrep(num2str(h),'.','p'),...
                'g',strrep(num2str(trained_net.g),'.','p'),...
                'seed',num2str(seed)]);
    
end
% hAx.CLim=[1 1.5]
%
% figure
% M=100;
% psi=(0:M-1)./M*2*pi;
% a=1;
% b=0.999;
% r=b*sqrt( cos(2*psi)+sqrt( (a/b)^4-sin(2*psi).^2   ) );
% z=[r.*cos(psi);r.*sin(psi)];
% plot(z(1,:),z(2,:),'o')