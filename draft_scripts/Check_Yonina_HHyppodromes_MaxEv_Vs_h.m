clear;
%close all;
ellipseFlag=1;
UniformFlag=0;
eigJW={};
rng(1);
if ellipseFlag
    h_vec=[-0.4 -0.2 0 0.2 0.4 0.8];  %  ellipse
    a=2;b=1;
    g_vec=[0.2];
else
    h_vec=[-0.3:0.1:0.2];  %  4hyppo
    a=4;b=1;
    g_vec=[1];
end

A=1.2;
% M_vec=[2:1:6 8 10 12 20 40 80]; % DUE TO TE BUG IT MUST BE DIVIDED BY 160
M_vec=[20]; % DUE TO TE BUG IT MUST BE DIVIDED BY 160

net.N=1000; net.phi=@(x)erf(x./sqrt(2)); net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.sim_resolution=240;

for gg=1:length(g_vec)
    gg
    tic;
    for jj=1:1
        net.g=g_vec(gg);
        [hp, net, sim] = prep_network_param(hp,net,struct);
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
%                 for pt=1:M_vec(mm)
                for pt=1

                    %                 ee_cl=eig((net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt))));
                    S=(net.W+net.wfb*net.wout')*diag(net.phip(x(:,pt)));
                    ee_cl=eigs(S,1,'largestreal');
                    eigJW{gg,hh,mm,pt,jj}=ee_cl;
                end
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
                
                lambda_cell{gg,hh,mm,jj}=lambda;
                
            end
        end
    end
    toc;
end

% TEST
% mm=1; pts=1+floor(hp.sim_resolution/M_vec(mm)/2)*[0:(M_vec(mm)-1)];plot(sim.f_ol(1,pts),sim.f_ol(2,pts),'o')
%% Load
IT=1;
for jj=1:IT
    jj
    for gg=1:length(g_vec)
        for hh=1:length(h_vec)
            hh
            for mm=1:length(M_vec)
%                 for pp=1:M_vec(mm)
for pp=1
                    %                 maxEvRSim(gg,hh,mm,pp,jj)=max(real(eigJW1{gg,hh,mm,pp,jj}-1));
                    maxEvPsiSim(gg,hh,mm,pp,jj)=max(real(eigJW{gg,hh,mm,pp,jj}-1));
                end
            end
            
        end
    end
end
%% Plot Max Lamda, Max over all points Vs M
colors=colormap;

for hh=1:length(h_vec)
    for mm=1:length(M_vec)
%         mLambda(hh,mm)=max(abs(squeeze(maxEvPsiSim(1,hh,mm,1:M_vec(mm)))));
        mLambda(hh,mm)=(abs(squeeze(maxEvPsiSim(1,hh,mm))));

    end
end

figure(11)
for hh=1:length(h_vec)
    semilogy(2*M_vec,mLambda(hh,:),'-o','Color',colors(hh*10,:));hold on;
    % plot(M_vec,mLambda(hh,:),'-o','Color',colors(hh*40,:));hold on;
    
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
xlim([4 40])
% title('Max of Max... abs Max EV Vs M- ellipse, g=0.2; Different M')

title('Max of Max... abs Max EV Vs M- 4Hyppo, g=1; Different M')
xlabel('2M'); ylabel('max Lambda (no abs)')
legend(leg)
box off

figure(12)
for hh=1:length(h_vec)
    semilogy(2*M_vec,mLambda(hh,:)./mLambda(hh,1),'-o','Color',colors(hh*10,:));hold on;
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
title(' Normalized abs Max EV Vs M- ellipse, g=0.2; Different M')
xlabel('2M'); ylabel('|max Lambda|')
legend(leg)
box off
%% Plot examples
colors=colormap;

for hh=1:length(h_vec)
    h=h_vec(hh)
sim.f_ol=A*1/(a-b)*[(a-b)*cos(sim.psi-rot)+h*cos((a-b)/b*(sim.psi-rot)); (a-b)*sin(sim.psi-rot)-h*sin((a-b)/b*(sim.psi-rot))]./(1+h/(a-b));

sim.f_ol=rotMat*sim.f_ol;

figure(11);plot(sim.f_ol(1,:),sim.f_ol(2,:),'-','Color',colors(hh*10,:)); hold on;
% hold on;plot(sim.f_ol(1,1),sim.f_ol(2,1),'x')
end
%% UP TO HERE!!

figure(1)
plot(h_vec,abs(squeeze(maxEvPsiSim(1,:,1,1:2))),'r-o');hold on;
plot(h_vec,abs(squeeze(maxEvPsiSim(1,:,2,1:4))),'g-o');
plot(h_vec,abs(squeeze(maxEvPsiSim(1,:,3,1:10))),'k-o')
title(' abs Max EV Vs h- 4hyppo, g=1; Different M')
xlabel('h'); ylabel('|max Lambda|')
legend('M=2(4)','M=4(8)','M=10(20)')
box off

figure(11)
semilogy(h_vec,abs(squeeze(maxEvPsiSim(1,:,1,1:2))),'r-o');hold on;
semilogy(h_vec,abs(squeeze(maxEvPsiSim(1,:,2,1:4))),'g-o');
semilogy(h_vec,abs(squeeze(maxEvPsiSim(1,:,3,1:10))),'k-o')
title('abs Max EV Vs h- 4hyppo, g=1, Log scale; Different M')
xlabel('h'); ylabel('log(|max Lambda|)')
legend('M=2(4)','M=4(8)','M=10(20)')
box off

figure(33)
% leg=num2str(h_vec(1));
title('Decay in correlation spectrum Vs m')
for hh=1:length(h_vec)
    leg{hh} = ['h = ',num2str(h_vec(hh))];
%     semilogy((lambda_cell{1,hh,4}),'-o')

    semilogy((lambda_cell{1,hh}),'-o')
hold on
end
legend(leg)
box off


%% Plot Vs M

figure(1)
for hh=1:length(h_vec)
    for mm=1:length(M_vec)
        semilogy(M_vec(mm),abs(squeeze(maxEvPsiSim(1,hh,mm,1:M_vec(mm)))),'o','Color',colors(hh*10,:));hold on;
    end
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
title(' abs Max EV Vs h- 4hyppo, g=1; Different M')
xlabel('M'); ylabel('|max Lambda|')
legend(leg)
box off
%%
colors=colormap;

for hh=1:length(h_vec)
    for mm=1:length(M_vec)
        mLambda(hh,mm)=max(abs(squeeze(maxEvPsiSim(1,hh,mm,1:M_vec(mm)))));
    end
end

figure(11)
for hh=1:length(h_vec)
    semilogy(2*M_vec,mLambda(hh,:),'-o','Color',colors(hh*10,:));hold on;
    % plot(M_vec,mLambda(hh,:),'-o','Color',colors(hh*40,:));hold on;
    
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
xlim([4 40])
% title('Max of Max... abs Max EV Vs M- ellipse, g=0.2; Different M')

title('Max of Max... abs Max EV Vs M- 4Hyppo, g=1; Different M')
xlabel('2M'); ylabel('max Lambda (no abs)')
legend(leg)
box off

figure(12)
for hh=1:length(h_vec)
    semilogy(2*M_vec,mLambda(hh,:)./mLambda(hh,1),'-o','Color',colors(hh*10,:));hold on;
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
title(' Normalized abs Max EV Vs M- ellipse, g=0.2; Different M')
xlabel('2M'); ylabel('|max Lambda|')
legend(leg)
box off




%%
M_vec=[ 2     3     4     5     6     8    10    12    20  ];
colors=colormap;
for hh=1:length(h_vec)
    for mm=1:length(M_vec)
        % mLambda(hh,mm)=mean(max((squeeze(maxEvPsiSim(1,hh,mm,1:M_vec(mm),:)))));
        mLambda(hh,mm)=mean(max(abs(squeeze(maxEvPsiSim(1,hh,mm,1:M_vec(mm),:)))));
        
    end
end

figure(11)
for hh=1:length(h_vec)
    semilogy(M_vec,mLambda(hh,:),'-o','Color',colors(hh*40,:));hold on;
    % plot(M_vec,mLambda(hh,:),'-o','Color',colors(hh*40,:));hold on;
    
    leg{hh} = ['h = ',num2str(h_vec(hh))];
end
% title('Max of Max... abs Max EV Vs M- ellipse, g=0.2; Different M')

title('Max of Max... abs Max EV Vs M- 4Hyppo, g=1; Different M')
xlabel('M'); ylabel('max Lambda (no abs)')
legend(leg)
box off
