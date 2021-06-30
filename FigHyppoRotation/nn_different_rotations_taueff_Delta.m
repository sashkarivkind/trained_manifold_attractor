% clear;
N=4000;
% g_vec=[1];%0:0.1:2;
g_vec=[1];%0:0.1:2;
% h_vec=[-0.2 0.2 0.4 0.6 0.8];
% h_vec=-0.15;
% h_vec=0.2;

% h_vec=0.2;
% DegRot=180;
% DegIC_vec=[0,90];
% FinalFlag=0; DegRot_vec=[182,180];


% % BAD 4Hyppo
h_vec=0.9;a=4;
% DegIC_vec=[0];
% FinalFlag=0; DegRot_vec=[90];
% DegIC_vec=[45];
% FinalFlag=0; DegRot_vec=[91];
% DegIC_vec=[20];
% FinalFlag=0; DegRot_vec=[160];
DegIC_vec=[0];
FinalFlag=0; DegRot_vec=[95];


% % 4Hyppo
% h_vec=-0.15;a=4;
% DegIC_vec=[315];
% FinalFlag=0; DegRot_vec=[136];
% DegIC_vec=[180];
% FinalFlag=0; DegRot_vec=[135];

% Ellipse
% h_vec=0.2;a=2;
% DegIC_vec=[0];
% FinalFlag=0; DegRot_vec=[182];
% DegIC_vec=[90];
% FinalFlag=0; DegRot_vec=[180];

% RING
% a=2;
% h_vec=0;
% DegRot=90;
% DegIC_vec=[270,300,0, 180,120];
% FinalFlag=1; DegRot_vec=[90 90 90 90 90];

b=1;
A=1.2;



% h_vec=[-0.5 0 0.2 0.4 0.6];
for seed=1%[1:4];%8:10;
    seed
    run_new_net = 1;
    zdelta_rec={};
    input_mode = 'hyppo';
    t_conv_vec=[];
    for gg=1:length(g_vec)
        for hh=1:length(h_vec)
            tic;
            gg
            if 1 %
                rng(seed);
                hp=struct;
                net=struct;
                hp.M=40;
                h=h_vec(hh);
                rot=0;
                sim=struct;
                
                net.g=g_vec(gg);
                net.N=N;
                
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
            end
            net.win=net.wfb;
            
            % Calculate tau_eff
            sim.x=x;
            tau_inv_vec=[];
            tau_eff2=[];
            
            all_pts=1:length(sim.f_ol);
            for ppp=1:length(all_pts)
                semi=semi_empirical_spectrum(net,sim,30,all_pts,ppp);
                % tau_inv_vec(end+1)=semi.tau_inv;
                tau_eff2(end+1)=semi.tau_approx_by_deriv;
                
            end
            % Calculate ev
 

            % Do rotations
            for ddiicc=1:length(DegIC_vec)
                DegIC=DegIC_vec(ddiicc);
                DegRot=DegRot_vec(ddiicc);
                [~,ptic]=min(abs(sim.psi-DegIC/180*pi));
                dt=1;
                tmax=2000;
                epsilonmag=0.01;
                tswitch=[50];
                pulse_width = [2000]/dt;
                if  FinalFlag
                    psi_vec=[DegRot]*pi/180;
                else
                    psi_vec=[DegIC+DegRot]*pi/180;
                end
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
%                 thet_ic=2*pi*rand(n_ic,1);
%                 psi_ic=[0:(net.N-1)]./net.N.*2*pi; 
                x_on_manifold_ic=x(:,ptic);
                
                tic;
                ptic
                [xdeltaRot,zdeltaRot{ddiicc},xrecdeltaRot]=fast_conv_to_fp_extended(net,[],...
                    struct('xinit',x_on_manifold_ic,...
                    'ol',0,...
                    'tmax',tmax,...
                    'dt',dt,...
                    'u_in',u_in,...
                    'save_neurons',[1:200:net.N]));
                toc;
            end
            % Calculate Delta
            x_epsi_on_manifold = fast_conv_to_fp(net,...
                sim.f_ol,...
                struct('ol_with_fixed_input',1,...
                'u',u_in(:,1000)*ones(1,size(sim.f_ol,2))));
            z_epsi_on_manifold = net.wout'*net.phi(x_epsi_on_manifold);
            psi_on_manifold = atan2(sim.f_ol(2,:),sim.f_ol(1,:));
            % plot(psi_on_manifold,atan2(z_epsi_on_manifold(2,:),z_epsi_on_manifold(1,:))-psi_on_manifold )
            Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
            
            %Save
                        save(['FigureResultsN',num2str(net.N),...
                            'a',strrep(num2str(a),'.','p'),...
                            'Deg',strrep(num2str(DegRot),'.','p'),...
                            'h',strrep(num2str(h),'.','p'),...
                            'g',strrep(num2str(net.g),'.','p'),...
                            'eps',strrep(num2str(epsilonmag),'.','p'),...
                            'seed',num2str(seed)]);
            toc;
        end
        
        
    end
end
%%
% figure;
% for ddiicc=1:length(DegIC_vec)
%     plot(zdeltaRot{ddiicc})
%     hold on;
% end
% plot(sim.f_ol(1,:),sim.f_ol(2,:),'--');
% 
% %%
% figure;
% for ddiicc=1:length(DegIC_vec)
%     %     plot(180/pi*unwrap(((angle(zdeltaRot{ddiicc})-angle(zdeltaRot{ddiicc}(1))))))
%     plot(unwrap(((angle(zdeltaRot{ddiicc})-angle(zdeltaRot{ddiicc}(1))))))
%     
%     %     plot((180/pi*(angle(zdeltaRot{ddiicc}))))
%     hold on;
% end