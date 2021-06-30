dt=1; itNum=1;
g_vec=[0.01 0.1:.1:0.9 0.95 1:0.1:1.8];
A=1.5;hp.A=A;
net.phi=@(x)erf(x./sqrt(2));
net.phip=@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.omega_vec=0;
results.G_polar={};
%
net.N=1000;
results.maxev=[];
hp.M=20;
% Calcualte the empirical max E.V. of the stability matrix in r directiom
for GG=1:length(g_vec)
    g=g_vec(GG);
    net.g=g;
    GG
    tic;
    for II=1:itNum
        result.with_learning={}; Gij={};
        hp.sim_resolution=hp.M;
        [hp, net, sim] = prep_network_param_learned_angles_only(hp, net,struct);
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        sim.r = net.phi(x);
        regfac = hp.M*hp.alpha_reg*eye(hp.M);
        net.wout = (sim.r/...
            (sim.r'*sim.r+regfac))*sim.f_ol'; %obtain least mean square solution
        cnt=0;
        pt=1;
        %         for pt = 1:hp.M
        cnt=cnt+1;
        results.evFull{GG,II}=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt)))-eye(net.N));

        results.ev{GG,II}=eig((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt)))-eye(net.N));
        %             results.bulk(GG,II)=max(real(eig((net.W)*diag(net.phip(x(:,pt)))-eye(net.N))));
        results.maxev(GG,II)=max(real(results.ev{GG,II}));
%         bulkNumeric(GG,II)=sqrt(net.g^2*(mean((net.phip(x(:,pt)).^2))))-1;
    end
    toc;
end
%% Calcualte the max E.V. of the bulk using mean field
mu_0=1;sig_0=1;
for gg=1:length(g_vec)
    g=g_vec(gg);
    [sigma]=solveFixedPointMuSigRing_Erf(net.phi,g,A,sig_0);
    bulkMF(gg)=bulk_fun(g,sigma,A,net.phip);
end
%%
% g_vec=[0.01 0.2:0.2:1.8];
% g_vec=[1:0.2:1.8];
% plot(g_vec(1:5),results.maxev(1:5,:),'ok ','MarkerSize',5);hold on;
% plot(g_vec(6:end),results.maxev(6:end,:),'og ','MarkerSize',5);hold on;

if itNum>1
   errorbar(g_vec,mean(results.maxev,2),std(results.maxev,[],2),'o','MarkerSize',5);hold on;
else
    plot(g_vec,results.maxev,'o','MarkerSize',5);hold on;
end
plot(g_vec,bulkMF,'-k')

% plot(g_vec(:),results.bulk(:,:),'or ','MarkerSize',5);hold on;
% plot(g_vec,mean(bulkNumeric,2),'--k')
% plot(g_vec,mean(bulkNumeric,2),'--k')

% Calcualte the max E.V. of in r direction using mean field

% MF lambda_r
g_vec=[0.01 0.1:.1:1.8];
ol_pole_vec=zeros(size(g_vec));
beta_info=load('beta_graph_long.mat');
for gg=1:length(g_vec)
    ii=find(abs(beta_info.J_vec-g_vec(gg))<1e-10);
    if length(ii)~=1
        error('pole info not found!');
    end
    lambdar(gg)=beta_info.debuu{ii}.spectrum_rr(1,1);
end

plot(g_vec,lambdar)
xlim([0 1.8])
ylim([-0.65 0.1])
axis square
box off
%% Plot the stability matrix
gg=5;
figure(20)
plot(results.evFull{gg,1},'.')
% xlim([-1.3 0.1]); ylim([-0.3 0.3])
xlim([-2.2 0.1]); ylim([-1 1])

title('g=0.4'); box off; axis square
gg=16;
figure(21)
plot(results.evFull{gg,1},'.')
xlim([-2.2 0.1]); ylim([-1 1])
title('g=1.4'); box off; axis square