clear;close all;
RunCLFlag=1;
UniformFlag=0;
beta_rescaled_rec={};
beta0_rescaled_rec={};
M_rescaled_rec={};
y_rescaled_rec={};
dd_rec={};

beta0_rec={};
M_rec={};
beta_rec={};
I28_debu={};
eigJW={};
% % h_vec=[-0.9:0.1:0.9];
h_vec=[0.15];

% h_vec=[ -0.2:0.04:0.2 ];
% h_vec=[0 0.02];
A=1.2;
%
% g_vec=[0.2 0.5 1  1.2 1.5];
M_vec=[40]; % DUE TO TE BUG IT MUST BE DIVIDED BY 160
% M_vec=[40]; % DUE TO TE BUG IT MUST BE DIVIDED BY 160

g_vec=[0.8];
% M_vec=[20 ];
net.N=1000; net.phi=@(x)erf(x./sqrt(2)); net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
a=4;b=1; K_vec=[24];
for gg=1:length(g_vec)
    gg
    tic;
    for jj=1:1
        net.g=g_vec(gg);
        [hp, net, sim] = prep_network_param(struct,net,struct);
        for mm=1:length(M_vec)
            mm
            hp.M=M_vec(mm); % BUG? DOES p gets updated in te relevant places??? FIX IT!!!
            
            for hh=1:length(h_vec)
                h=h_vec(hh);
                rot=0;
                rotMat=[cos(rot), -sin(rot); sin(rot), cos(rot) ];
                if UniformFlag
                    tmp_f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                    sim.f_ol=curvspace(tmp_f_ol',hp.sim_resolution);
                    sim.f_ol=sim.f_ol';
                else
                    sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                end
                sim.f_ol=rotMat*sim.f_ol;
                
                figure(11);plot(sim.f_ol(1,:),sim.f_ol(2,:),'o'); hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
                
                x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1)); % Sould update before the resolution to be hp.M!!!
                
                %learning output weights
                sim.r = net.phi(x);
                pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
                regfac = hp.M*hp.alpha_reg*eye(length(pts));
                net.wout = (sim.r(:,pts)/...
                    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
                
                %obtaining open loop output
                sim.z_ol = net.wout'*sim.r;
                %obtaining the eigenvalues
                pt=1;
                ee_cl=eigs((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))),1,'largestreal');
                ee_cl1=eigs((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt))),1,'largestreal');
                ee_cl2=eigs((net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(x(:,pt))),1,'largestreal');
                
                
                    
                    
                eigJW{gg,hh,mm,jj}=ee_cl; eigJW1{gg,hh,mm,jj}=ee_cl1; eigJW2{gg,hh,mm,jj}=ee_cl2;
                
                %obtaining the semi-MF eigenvalues
                pts=1:1:length(sim.f_ol);
                Meff=length(pts);
                r=sim.r(:,pts);
                
                g=net.g;  sim.rp=net.phip(x); rp=sim.rp(:,pts);
                x1full=(x-net.wfb*sim.f_ol);  x1=x1full(:,pts);
                
                Crr=r'*r;  [eta,lambda]=eig(Crr);
                lambda=1/Meff*1/net.N*lambda; %doublecheck normalization
                
                [lambda,ii]=sort(diag(lambda),'descend');
                lambda=abs(real(lambda))+1e-20;
                eta=eta(:,ii); eta0=diag(eta(:,1));
                
                for kk=1:length(K_vec)
                    K=K_vec(kk);
                    lambdaK=lambda(1:K);  etaK=eta(:,1:K);
                    
                    an=x1*etaK*inv(diag(net.g*sqrt(Meff*lambdaK))); %doublecheck normalization
                    zzz=1;
                    
                    % equation
                    for dd=1:2
                        tmp=(g*diag(sqrt(Meff*lambdaK)));
                        lhs=tmp^2;
                        B0=g^2*1/net.N*(diag(rp(:,pt))*(net.wfb(:,dd)))'*(r*etaK);
                        B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*etaK);
                        B0=B0';
                        B1=B1'*tmp;
                        alpha= inv(zzz*lhs-B1)*B0/zzz;
                        % recovering gain
                        rpc=r*etaK;
                        qn=net.wout'*rpc/diag(lambdaK)/Meff;
                        % qn_alt=net.wout'*rpc(rpc*rpc');
                        Gnprefac=(g^(-1)*sqrt(Meff*lambdaK')*tmp)';
                        Gn=zzz*Gnprefac.*alpha;
                        MM=inv(lhs)* (B1+B0*((qn(dd,:).*Gnprefac')));
                        %             MM_m=(B1+B0*((qn(dd,:).*Gnprefac')));
                        G=qn*Gn;
                        y=eig(MM)-1;
                        cp=(rp'*rp);
                        I28_debu{gg,hh,mm,jj,dd,kk}=etaK'*diag(cp(1,:))*etaK*diag(lambdaK);
                        beta0_rescaled_rec{gg,hh,mm,jj,dd,kk}=inv(lhs)*B0*((qn(dd,:).*Gnprefac'));
                        M_rescaled_rec{gg,hh,mm,jj,dd,kk}=MM;
                        y_rescaled_rec{gg,hh,mm,jj,dd,kk}=y;
                        beta_rescaled_rec{gg,hh,mm,jj,dd,kk}=inv(lhs)*B1;
                        dd_rec{gg,hh,mm,jj,dd}=dd;
                    end
                end
                
                if RunCLFlag
                    %simulating closed loop
                    x_test_rand=fast_conv_to_fp(net,[],struct('xinit',5*randn(size(x))));
                    r_test_rand=net.phi(x_test_rand);
                    sim.z_test_rand= net.wout'*r_test_rand;
                    
                    x_test_noi=fast_conv_to_fp(net,[],struct('xinit',x+1e-3*randn(size(x))));
                    r_test_noi=net.phi(x_test_noi);
                    sim.z_test_noi= net.wout'*r_test_noi;
                    
                    %plotting results
                    figure(1000+gg); 
%                     subplot(5,5,hh)
                    plot(sim.z_ol(1,:),sim.z_ol(2,:),'.');
                    hold on;
                    plot(sim.f_ol(1,:),sim.f_ol(2,:),'o');
                    plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'b+','linewidth',3);
                    plot(sim.z_test_rand(1,:),sim.z_test_rand(2,:),'x', 'linewidth',1);
                    plot(sim.z_test_noi(1,:),sim.z_test_noi(2,:),'d', 'linewidth',1);
                    
                    xlim([-2 2]); ylim([-2 2])
                end
            end
        end
    end
    toc;
end
%% Load
IT=1;
for jj=1:IT
    jj
    for gg=1:length(g_vec)
        for hh=1:length(h_vec)
            for mm=1:length(M_vec)
                maxEvRSim(gg,hh,mm,jj)=max(real(eigJW1{gg,hh,mm,jj}-1));
                maxEvPsiSim(gg,hh,mm,jj)=max(real(eigJW2{gg,hh,mm,jj}-1));
                for kk=1:length(K_vec)
                    maxEvR(gg,hh,mm,jj,kk)=max(real(y_rescaled_rec{gg,hh,mm,jj,1,kk}));
                    maxEvPsi(gg,hh,mm,jj,kk)=max(real(y_rescaled_rec{gg,hh,mm,jj,2,kk}));
                end
            end
            
        end
    end
end
%% Plot Lambdapsi Vs M, diff h M_vec=[2 4  10 16 20]; g=1.2
IT=1;     gg=3; g_vec(gg)
colors=colormap;
figure(50+gg)
for hh=1:length(h_vec)
    for jj=1:IT
        semilogy(M_vec,squeeze(abs(maxEvPsiSim(gg,hh,:))),'-o','Color',colors(2*hh,:));hold on
    end
end
title(['g=' num2str(g_vec(gg))])
box off;xlabel('M'); ylabel('MaxReal \lambda in R')

% plot the MF
% for jj=1:IT
%     for hh=1:length(h_vec)
%         for kk=[ 9]%:length(K_vec)]
%             semilogy(M_vec,squeeze(abs(maxEvPsi(gg,hh,:,jj,kk))),'-','Color',colors(hh,:)); hold on
%         end
%         %             ylim([-0.65 2])
%         % title(['h=' num2str(h_vec(hh))])
%     end
% end
    %% Plot Lambdapsi Vs h, diff M M_vec=[2 4  10 16 20]; g=1.2
    dd=2;logFlag=0;
    gg=1;    colors=colormap;
    kmax= find(K_vec==2);
    figure(60+gg)
    %     for gg=1:length(g_vec)
    %         subplot(3,2,gg)
    for mm=length(M_vec)
        for jj=1:1
            if dd==1
                plot(h_vec,squeeze(maxEvRSim(gg,:,mm,jj)),'o','Color',colors(mm*10,:))
            else
                if logFlag
                    semilogy(h_vec,squeeze(abs(maxEvPsiSim(gg,:,mm,jj))),'o','Color',colors(mm*10,:))
                else
                    plot(h_vec,squeeze((maxEvPsiSim(gg,:,mm,jj))),'o','Color',colors(mm*10,:))
                    
                end
            end
            
            hold on
        end
    end
    title(['g=' num2str(g_vec(gg))])
    box off
    xlabel('h'); ylabel('MaxReal \lambda in R')
    %     end
    %                 legend('M=4 (8)',     '8'    ,'10'    ,'16'    ,'20')
    % %
    % plot the MF
    %     for gg=1:length(g_vec)
    %         subplot(3,2,gg)
    %%
    for jj=1:1
        for mm=length(M_vec)
            ind=0;
            for kk=1:2:length(K_vec)
                ind=ind+1;
              leg{ind} = ['K = ',num2str(K_vec(kk))];

                if dd==1
                    plot(h_vec,squeeze(maxEvR(gg,:,mm,jj,kk)),'--','Color',colors(kk*5,:)); hold on
                else
                    if logFlag
                    semilogy(h_vec,squeeze(abs(maxEvPsi(gg,:,mm,jj,kk))),'--','Color',colors(kk*5,:)); hold on
                    else
                     plot(h_vec,squeeze((maxEvPsi(gg,:,mm,jj,kk))),'--','Color',colors(kk*5,:)); hold on
                    end
                end
                
            end
            %             ylim([0 .05])
        end
        %         end
%         title(['g=' num2str(g_vec(gg)) '   K=' num2str(K_vec(kk))])
    end
    %     xlim([-0.5 0.5])
    xlim([-0.9 0.9])
    legend(leg)
%     legend('M=4 (8)','M=10 (20)','M=16 (32)','M=20 (40)')
%         legend('M=4 (8)','M=10 (20)','M=16 (32)','M=20 (40)','M=40 (80)')

%         legend('M=2 (4)','M=4 (8)','M=10 (20)','M=16 (32)','M=20 (40)','M=40 (80)')

% xlim([-0.8 0.9])
% 
% xlim([-0.2 0.2])
% ylim([-0.025 0.025])
    %% g HELPS
      %% Plot Lambdapsi Vs h, diff M M_vec=[2 4  10 16 20]; g=1.2
    dd=2;logFlag=1;
    gg=2;    colors=colormap;

    figure(60+gg)
    %     for gg=1:length(g_vec)
    %         subplot(3,2,gg)
    for gg=1:5%length(g_vec)
    for mm=6
        for jj=1:1
            if dd==1
                plot(h_vec,squeeze(maxEvRSim(gg,:,mm,jj)),'-o','Color',colors(gg*10,:))
            else
                if logFlag
                    semilogy(h_vec,squeeze(abs(maxEvPsiSim(gg,:,mm,jj))),'-o','Color',colors(gg*10,:))
                else
                    plot(h_vec,squeeze((maxEvPsiSim(gg,:,mm,jj))),'-o','Color',colors(gg*10,:))
                    
                end
            end
            
            hold on
        end
    end
    end
    plot(h_vec,zeros(size(h_vec)),'-k')
%     title(['g=' num2str(g_vec(gg))])
legend('g=0.5','g=1','g=1.2','g=1.5')
    box off
    xlabel('h'); ylabel('MaxReal \lambda in R')
    %%
        set(gca,'ColorOrderIndex',1);

    kmax= find(K_vec==4);
    for gg=1:2%length(g_vec)
    for jj=1:1
        for mm=5
            for kk=[kmax]%:length(K_vec)]
                if dd==1
                    plot(h_vec,squeeze(maxEvR(gg,:,mm,jj,kk)),'--','Color',colors(gg*10,:)); hold on
                else
                    if logFlag
                    semilogy(h_vec,squeeze(abs(maxEvPsi(gg,:,mm,jj,kk))),'--','Color',colors(gg*10,:)); hold on
                    else
                     plot(h_vec,squeeze((maxEvPsi(gg,:,mm,jj,kk))),'--','Color',colors(gg*10,:)); hold on
                    end
                end
                
            end
            %             ylim([0 .05])
        end
        %         end
%         title(['g=' num2str(g_vec(gg)) '   K=' num2str(K_vec(kk))])
    end
    end
    %     xlim([-0.5 0.5])
        plot(h_vec,zeros(size(h_vec)),'-k')
%     title(['g=' num2str(g_vec(gg))])
legend('g=0.2','g=0.5','g=1','g=1.2','g=1.5')
    xlim([-0.9 0.9])
