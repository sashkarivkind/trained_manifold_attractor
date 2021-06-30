clear;close all;
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
h_vec=[0.1 0.6 0.7 0.9];
for hh=1:4
    hh
    a=6;b=1;h=h_vec(hh);
    K=24;
    net.g=0.6;
    rot=pi/a;
    % rot=pi/4;
    
    % rot=0;
    net.N=1000;
    
    for jj=1:1
        
        jj
        % dd=2; %direction - 1 - transversal, 2- parallel
        
        net.phi=@(x)erf(x./sqrt(2));
        net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
        [hp, net, sim] = prep_network_param(struct,net,struct);
        % sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];
        % rot=pi/a;
        
        rotMat=[cos(rot), -sin(rot); sin(rot), cos(rot) ];
        sim.f_ol=1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))];
        
        % % % A=1.2;
        % % % f_ol_circ=A*[cos(sim.psi-rot);sin(sim.psi-rot)];
        % % % f_ol_rhomb=A*[abs(cos(sim.psi-rot)).*cos(sim.psi-rot);abs(sin(sim.psi-rot)).*sin(sim.psi-rot)];
        % % %                 sim.f_ol=(1-h)*f_ol_circ+h*f_ol_rhomb; % b now measures how much square is introduced
        
        
        sim.f_ol=rotMat*sim.f_ol;figure(11);plot(sim.f_ol(1,:),sim.f_ol(2,:),'o'); hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
        hp.M=20;
        
        x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));
        
        %learning output weights
        sim.r = net.phi(x);
        pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
        regfac = hp.M*hp.alpha_reg*eye(length(pts));
        net.wout = (sim.r(:,pts)/...
            (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
        
        %obtaining open loop output
        sim.z_ol = net.wout'*sim.r;
        
        
        
        
        
        
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
        
        
        
        
        
        
        
        pt=1;
        
        
        
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
            
            MM_m=(B1+B0*((qn(dd,:).*Gnprefac')));
            G=qn*Gn
            y=eig(MM)-1;
            
            cp=(rp'*rp);
            I28_debu{hh,jj,dd}=eta'*diag(cp(1,:))*eta*diag(lambda);
            
            beta0_rescaled_rec{hh,jj,dd}=inv(lhs)*B0*((qn(dd,:).*Gnprefac'));
            M_rescaled_rec{hh,jj,dd}=MM;
            y_rescaled_rec{hh,jj,dd}=y;
            ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
            eigJW{hh,jj,dd}=ee_cl;
            
            beta_rescaled_rec{hh,jj,dd}=inv(lhs)*B1
            
            
            % beta0Only_rec{end+1}=B0;
            
            
            beta0_rec{hh,jj,dd}=B0*((qn(dd,:).*Gnprefac'));
            M_rec{hh,jj,dd}=MM_m;
            beta_rec{hh,jj,dd}=B1;
            dd_rec{hh,jj,dd}=dd;
            
        end
    end
end
%%%%
% lhs=1;
%%
IT=1;
close all
% vecind=2:2:K;
% vecind=1:2:K;
% vecind=2:2:6;
vecind=2:2:6;
b0=[];b1=[];Mmat=[];MmatNull=[];b0null=[];
dd=1 ;
for hh=1:4
    for jj=1:IT
%             b0(hh,jj,:,:)=log10(abs(beta0_rescaled_rec{hh,jj,dd}));
%             b1(hh,jj,:,:)=log10(abs(beta_rescaled_rec{hh,jj,dd}));
%             Mmat(hh,jj,:,:)=log10(abs(M_rescaled_rec{hh,jj,dd}));
            b0(hh,jj,:,:)=(beta0_rescaled_rec{hh,jj,dd});
            b1(hh,jj,:,:)=(beta_rescaled_rec{hh,jj,dd});
            Mmat(hh,jj,:,:)=(M_rescaled_rec{hh,jj,dd});
            
            b0null(hh,jj,:,:)=beta0_rescaled_rec{hh,jj,dd};
            b0null(hh,jj,:,3:end)=zeros(size(b0null(hh,jj,:,3:end)));
            MmatNull(hh,jj,:,:)=b1(hh,jj,:,:)+b0null(hh,jj,:,:);
    end
    figure(10) 
    subplot (4,2,2*hh-1)
    plot(eigJW{hh,jj,1}-1,'*');hold on;
        xlim([-1.4 0.1])

    subplot (4,2,2*hh)
    plot(y_rescaled_rec{hh,jj,dd},'x','MarkerSize',15);hold on;
    plot(y_rescaled_rec{hh,jj,dd},'xr','MarkerSize',15);hold on;
    legend('trans','tangent')
    xlim([-1.4 0.1])
    title(num2str(h_vec(hh)))
    
    figure(20)
    subplot (4,1,hh)
    imagesc(squeeze(mean(b0(hh,:,vecind,vecind),2)));
    colorbar
    title(['mean beta0,' num2str(h_vec(hh))])
    
    figure(21)
    subplot (4,1,hh)
    imagesc(squeeze(mean(b1(hh,:,vecind,vecind),2)));
    colorbar
    title(['mean beta1,' num2str(h_vec(hh))])

    figure(22)
    subplot (4,1,hh)
    imagesc(squeeze(mean(Mmat(hh,:,vecind,vecind),2)));
    colorbar
    title(['mean M,' num2str(h_vec(hh))])
       
    figure(50)
    subplot (4,1,hh)
    for jj=1:IT
        mtmp=M_rescaled_rec{hh,jj,dd};
        [vv,ee]=eig(mtmp(vecind,vecind)-eye(length(vecind)),'vector');
        plot(real(vv(:,1)));hold on
        ee(1)
    title(['PC1,' num2str(h_vec(hh))])
    end
    
        figure(51)
    subplot (4,1,hh)
    for jj=1:IT
        mtmp=squeeze(MmatNull(hh,jj,:,:));
        [vv2,ee2]=eig(mtmp(vecind,vecind)-eye(length(vecind)),'vector');
        plot(real(vv2(:,1)));hold on
        ee2(1)
    title(['PC1,' num2str(h_vec(hh))])
    end
end

%%
for qq=1:4
    figure(201)
    
    subplot (4,1,qq)
    plot(eigJW{2*qq-1}-1,'*');hold on;
    % xlim([-1.4 0.1])
    figure(301)
    subplot (4,1,qq)
    plot(y_rescaled_rec{2*qq-1},'x','MarkerSize',15);hold on;
    plot(y_rescaled_rec{2*qq},'xr','MarkerSize',15);hold on;
    legend('trans','tangent')
    % xlim([-1.4 0.1])
    title(num2str(h_vec(qq)))
end
%%
figure;
for qq=1:4
    figure(21)
    subplot (4,1,qq)
    
    imagesc(log10(abs(beta0_rescaled_rec{2*qq-1})));%,[-2 2]);
    colorbar
    title(num2str(h_vec(qq)))
    figure(22)
    subplot (4,1,qq)
    
    imagesc(log10(abs(beta_rescaled_rec{2*qq-1})));%,[-2 2]);
    title(num2str(h_vec(qq)))
    
    colorbar
    figure(23)
    subplot (4,1,qq)
    imagesc(log10(abs(M_rescaled_rec{2*qq-1})));%,[-2 2]);
    colorbar
    title(num2str(h_vec(qq)))
    
    
end
%%
figure;
plot(ee_cl-1,'*');
hold on;hold on;plot(y,'x','MarkerSize',15)
figure;imagesc(log10(abs(MM_m)))
figure;imagesc(log10(abs(inv(lhs)*B1)));title('beta1')
figure;imagesc(log10(abs(inv(lhs)*B0*((qn(dd,:).*Gnprefac')))));title('beta0')
%%
zz=[];
for qq=1:1
    zz(qq,:,:)=log10(abs(beta0_rescaled_rec{qq}));
    %     zz(qq,:,:)=(abs(I28_debu{qq}));
    
    %     zz(qq,:,:)=log10(abs(beta_rec{qq}));
end
% mean_mat=reshape(mean(zz,1),K,K);
% std_mat=reshape(std(zz,[],1),K,K)
mean_mat=squeeze(mean(zz,1));
std_mat=squeeze(std(zz,[],1));

figure;
subplot(2,2,1);
imagesc(mean_mat)
title('mean')
colorbar
subplot(2,2,2);
imagesc(std_mat)
title('std')
colorbar

subplot(2,2,3);
imagesc(mean_mat-diag(diag(mean_mat)))
title('mean no diag')
colorbar

subplot(2,2,4);
imagesc(std_mat-diag(diag(std_mat)))
title('std no diag')
colorbar

