clear;
grid_res=0.025;
N=4000;
% g_vec=[1];%0:0.1:2;
g_vec=[1];%0:0.1:2;
h_vec=[0];
A=1.2;
a=4;b=1;
DegRot=30;
% h_vec=[-0.5 0 0.2 0.4 0.6];
for seed=[1:3];%8:10;
    seed
    run_new_net = 1;
    zdelta_rec={};
    % g_vec=0:0.2:2;
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
            
            dt=0.1;
            tmax=3000;
            epsilonmag=0.01;
            tswitch=[50];
            pulse_width = [2000]/dt;
            psi_vec=[DegRot]*pi/180;
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
            [xdelta,zdelta,xrecdelta]=fast_conv_to_fp_extended(net,[],...
                struct('xinit',x(:,1),...
                'ol',0,...
                'tmax',tmax,...
                'dt',dt,...
                'u_in',u_in,...
                'save_neurons',[1:200:net.N]));
            %%
            %             psi_hat=psi_decoded([real(zdelta);imag(zdelta)],sim.f_ol,sim.psi);
            %             figure;
            %             subplot(2,1,1);
            %             plot(angle(zdelta));
            %             hold on;
            %             plot(psi_hat);
            %             hold on;
            %             plot(abs(zdelta))
            %             subplot(2,1,2);
            %             plot(xrecdelta');
            %             %%
            %             figure;
            %             subplot(2,1,1);
            %             plot(180/pi*psi_hat);
            %             hold on;
            %             plot(abs(zdelta))
            %             subplot(2,1,2);
            %             plot(xrecdelta');
            %             %%
            q=0.01;
            t_conv=abs(angle(zdelta)-angle(zdelta(end)))<q*abs(angle(zdelta(1))-angle(zdelta(end)));
            t_conv_vec(end+1)=min(find(t_conv))*dt-tswitch;
            zdelta_rec{gg}=zdelta;
            
            z_vec=-2:grid_res:2;
            rr=1;
            z_grid = repmat(z_vec,length(z_vec),1);
            z1=z_grid(:);
            z_grid_t=z_grid';
            z2=z_grid_t(:);
            f_ol_tag=[z1';z2'];
            psi_aux=atan2(f_ol_tag(2,:),f_ol_tag(1,:));
            rho_aux=abs(f_ol_tag(2,:)+1i*f_ol_tag(1,:));
            
            ext_input= u_in(:,round(tswitch/dt));
%             tic;
            x_epsi0 = fast_conv_to_fp(net,...
                f_ol_tag,...
                struct('ol_with_fixed_input',1,...
                'u',0*ext_input*ones(1,size(f_ol_tag,2))));
            z_epsi0 = net.wout'*net.phi(x_epsi0);
            
            x_epsi = fast_conv_to_fp(net,...
                f_ol_tag,...
                struct('ol_with_fixed_input',1,...
                'u',ext_input*ones(1,size(f_ol_tag,2))));
            z_epsi = net.wout'*net.phi(x_epsi);
%             toc;
            
            % net.W=[];
            % net.wout=[];
            % net.win=[];

             save(['2DPlotsN',num2str(net.N),...
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