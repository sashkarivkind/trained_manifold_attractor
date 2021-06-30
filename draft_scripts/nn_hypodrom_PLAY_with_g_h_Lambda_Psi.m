clear;close all;
RunCLFlag=0;
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
% h_vec=[-0.9:0.1:0.9];
h_vec=[-0.9:0.1:-0.3 -0.2:0.04:0.2 0.3:0.1:0.9];
% h_vec=[0 0.02];
A=1.5;
%
% g_vec=[0.2 0.5 1  1.2 1.5];
% M_vec=[4 8 10 16 20 ];
g_vec=[1];
M_vec=[20 ];
% M_vec=[20];

%
for gg=1:length(g_vec)
    gg
    for jj=1:5
        net.N=1000;    net.g=g_vec(gg);
        
        net.phi=@(x)erf(x./sqrt(2));
        net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
        [hp, net, sim] = prep_network_param(struct,net,struct);
        for mm=1:length(M_vec)
            mm
            hp.M=M_vec(mm);
            for hh=1:length(h_vec)
                %                 hh
                a=6;b=1;h=h_vec(hh);
                K=24;
                
                % rot=pi/4;
                
                rot=0;
                
                % dd=2; %direction - 1 - transversal, 2- parallel
                
                % sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];
                % rot=pi/a;
                
                rotMat=[cos(rot), -sin(rot); sin(rot), cos(rot) ];
                %         sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];

                if UniformFlag
                    tmp_f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                    sim.f_ol=curvspace(tmp_f_ol',hp.sim_resolution);
                    sim.f_ol=sim.f_ol';
                    % % % A=1.2;
                    % % % f_ol_circ=A*[cos(sim.psi-rot);sin(sim.psi-rot)];
                    % % % f_ol_rhomb=A*[abs(cos(sim.psi-rot)).*cos(sim.psi-rot);abs(sin(sim.psi-rot)).*sin(sim.psi-rot)];
                    % % %                 sim.f_ol=(1-h)*f_ol_circ+h*f_ol_rhomb; % b now measures how much square is introduced
                else
                    sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                end
                sim.f_ol=rotMat*sim.f_ol;
                
                figure(11);plot(sim.f_ol(1,:),sim.f_ol(2,:),'o'); hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
                %         =20;
                
                x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
                
                %learning output weights
                sim.r = net.phi(x);
                pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
                regfac = hp.M*hp.alpha_reg*eye(length(pts));
                net.wout = (sim.r(:,pts)/...
                    (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
                
                %obtaining open loop output
                sim.z_ol = net.wout'*sim.r;
                
                
                pt=1;
                ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
                ee_cl1=eig((net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(x(:,pt))));
                ee_cl2=eig((net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(x(:,pt))));

                eigJW{gg,hh,mm,jj}=ee_cl;
                eigJW1{gg,hh,mm,jj}=ee_cl1;
                eigJW2{gg,hh,mm,jj}=ee_cl2;
                
                pts=1:1:length(sim.f_ol);
                Meff=length(pts);
                r=sim.r(:,pts);
                
                g=net.g;
                sim.rp=net.phip(x);
                rp=sim.rp(:,pts);
                x1full=(x-net.wfb*sim.f_ol);
                x1=x1full(:,pts);
                
                Crr=r'*r;
                [eta,lambda]=eig(Crr);
                lambda=1/Meff*1/net.N*lambda; %doublecheck normalization
                
                [lambda,ii]=sort(diag(lambda),'descend');
                lambda=abs(real(lambda))+1e-20;
                eta=eta(:,ii);
                eta0=diag(eta(:,1));
                
                lambda=lambda(1:K);
                eta=eta(:,1:K);
                
                an=x1*eta*inv(diag(net.g*sqrt(Meff*lambda))); %doublecheck normalization
                zzz=1;
                
                % equation
                for dd=1:2
                    
                    tmp=(g*diag(sqrt(Meff*lambda)));
                    lhs=tmp^2;
                    B0=g^2*1/net.N*(diag(rp(:,pt))*(net.wfb(:,dd)))'*(r*eta);
                    B1=g^2*1/net.N*(diag(rp(:,pt))*(an))'*(r*eta);
                    B0=B0';
                    B1=B1'*tmp;
                    alpha= inv(zzz*lhs-B1)*B0/zzz;
                    % recovering gain
                    rpc=r*eta;
                    qn=net.wout'*rpc/diag(lambda)/Meff;
                    % qn_alt=net.wout'*rpc(rpc*rpc');
                    Gnprefac=(g^(-1)*sqrt(Meff*lambda')*tmp)';
                    Gn=zzz*Gnprefac.*alpha;
                    MM=inv(lhs)* (B1+B0*((qn(dd,:).*Gnprefac')));
                    %             MM_m=(B1+B0*((qn(dd,:).*Gnprefac')));
                    G=qn*Gn
                    y=eig(MM)-1;
                    cp=(rp'*rp);
                    I28_debu{gg,hh,mm,jj,dd}=eta'*diag(cp(1,:))*eta*diag(lambda);
                    beta0_rescaled_rec{gg,hh,mm,jj,dd}=inv(lhs)*B0*((qn(dd,:).*Gnprefac'));
                    M_rescaled_rec{gg,hh,mm,jj,dd}=MM;
                    y_rescaled_rec{gg,hh,mm,jj,dd}=y;
                    beta_rescaled_rec{gg,hh,mm,jj,dd}=inv(lhs)*B1
                    %
                    %             beta0_rec{hh,jj,dd}=B0*((qn(dd,:).*Gnprefac'));
                    %             M_rec{hh,jj,dd}=MM_m;
                    %             beta_rec{hh,jj,dd}=B1;
                    dd_rec{gg,hh,mm,jj,dd}=dd;
                    
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
                    figure(1000+gg); subplot(5,5,hh)
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
end
    %%
for    jj=1:5
    jj
    for gg=1:length(g_vec)
        for hh=1:length(h_vec)
            for mm=1:length(M_vec)
% % % % % % %                 figure(gg)
% % % % % % %                 
% % % % % % %                 subplot (5,5,mm)
% % % % % % %                 plot(eigJW{gg,hh,mm,jj}-1,'*');hold on;
% % % % % % %                 
% % % % % % %                 %     xlim([-1.4 0.2])
% % % % % % %                 title(num2str(h_vec(hh)))
                %     figure(301)
                maxEvRSim(gg,hh,mm,jj)=max(real(eigJW2{gg,hh,mm,jj}-1));
                
                maxEvR(gg,hh,mm,jj)=max(real(y_rescaled_rec{gg,hh,mm,jj,2}));
                %     subplot (5,5,hh)
% % % % % % % % % % % %                 plot(y_rescaled_rec{gg,hh,mm,jj,2},'xr','MarkerSize',15);hold on;
% % % % % % % % % % % %                 
% % % % % % % % % % % %                 plot(y_rescaled_rec{gg,hh,mm,jj,1},'xg','MarkerSize',15);hold on;
% % % % % % % % % % % %                 %     legend('','trans','tangent')
% % % % % % % % % % % %                     xlim([-0.05 0.05])
% % % % % % % % % % % %                     ylim([-0.05 0.05])
% % % % % % % % % % % %                 title(num2str(h_vec(hh)))
                %
                %     figure(401)
                %     subplot (4,3,qq)
                %     h=h_vec(qq);
                %             sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));
                %         sim.f_ol=rotMat*sim.f_ol;plot(sim.f_ol(1,:),sim.f_ol(2,:),'o'); hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
                %     xlim([-2 2]); ylim ([-2 2])
                %         title(num2str(h_vec(qq)))
            end
            
        end
    end
end
    %%
    colors=colormap;
    figure(55)
    for gg=1:length(g_vec)
        subplot(3,2,gg)
        for mm=1:length(M_vec)
            for jj=1:5
            plot(h_vec,squeeze(maxEvRSim(gg,:,mm,jj)),'o','Color',colors(mm*10,:))
            hold on
            end
        end
        title(['g=' num2str(g_vec(gg))])
        box off
        xlabel('h'); ylabel('MaxReal \lambda in \psi')
    end
                legend('M=4 (8)',     '8'    ,'10'    ,'16'    ,'20')
%%
    % plot the MF
    for gg=1:length(g_vec)
        subplot(3,2,gg)
        for jj=1:5
        for mm=1:length(M_vec)
            plot(h_vec,squeeze(maxEvR(gg,:,mm,jj)),'--','Color',colors(mm*10,:)); hold on
        end
        end
        title(['g=' num2str(g_vec(gg))])
    end
    
    %% Plot lambda Vs M for different h
    colors=colormap;
    figure(55)
    for gg=1:length(g_vec)
        subplot(3,2,gg)
        for hh=[1 5 13 20 25]
            for jj=1:1
            plot(M_vec,squeeze(maxEvRSim(gg,hh,:,jj)),'-o','Color',colors(hh*2,:))
            hold on
            end
        end
        ylim([-0.15 0.15])
        title(['g=' num2str(g_vec(gg))])
        box off
        xlabel('M'); ylabel('MaxReal \lambda in \psi')
    end
                legend('h=-0.9',     '-0.5'    ,'0'    ,'0.4'    ,'0.9')
                