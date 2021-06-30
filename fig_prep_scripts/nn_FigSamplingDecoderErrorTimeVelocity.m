clear;close all
saveFlag=0;
result=struct;
tmax=1000;
dt=1;
rng(1)

net.g=1.5;

net.N=1000;
net.phi=@(x)erf(x./sqrt(2));
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);
hp.sim_resolution=120;

[hp, net, sim] = prep_network_param(hp, net,struct);
% net = update_net_g(net,net.g); % this is a source for bugs and you generate two random matrices.
%
x = fast_conv_to_fp(net,sim.f_ol,struct('ol',1));

sim.r = net.phi(x);
size(sim.r )
result.with_learning={};
M_vec=[3,6,20]; % sim_resolution must be a multiplier of 2M
for MM=1:length(M_vec)
    hp.M=M_vec(MM);
    
    pts=1+floor(hp.sim_resolution/hp.M/2)*[0:(hp.M-1)]; % for full circle
    regfac = hp.M*hp.alpha_reg*eye(length(pts));
    net.wout = (sim.r(:,pts)/...
        (sim.r(:,pts)'*sim.r(:,pts)+regfac))*sim.f_ol(:,pts)'; %obtain least mean square solution
    
    result.with_learning{MM}.z_ol = net.wout'*sim.r;
    [x_cl,result.with_learning{MM}.z_cl]=fast_conv_to_fp_extended(net,[],...
        struct('xinit',x,...
        'ol',0,...
        'tmax',tmax,...
        'dt',dt));
    
   
    
end

%%
for MM=1:3
vel{MM}=mean(diff(mod(angle(result.with_learning{MM}.z_cl(:,1:2)),2*pi),[],2),2);
end

Er=phase([1,1i]*result.with_learning{1}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on; 
plot(sim.psi-pi,vel{1},'--r');
M=M_vec(1)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')


box off;
%ylim([-0.03 0.03]);xlim([-pi pi])



if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM6g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM6g0_5rnd1N2k.pdf')
    end
end

Er=phase([1,1i]*result.with_learning{2}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on; 
plot(sim.psi-pi,vel{2},'--r');
M=M_vec(2)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')
box off;%ylim([-0.03 0.03]);
xlim([-pi pi])


if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM12g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM12g0_5rnd1N2k.pdf')
    end
end

Er=phase([1,1i]*result.with_learning{3}.z_ol)-sim.psi;
figure;plot(sim.psi-pi,Er);
hold on; 
plot(sim.psi-pi,vel{3},'--r');
M=M_vec(3)*2;pts=1+floor(hp.sim_resolution/M)*[0:(M-1)];
plot(sim.psi-pi,zeros(size(sim.psi)),'k')
plot(sim.psi(pts)-pi,Er(pts),'.')
box off;%ylim([-0.03 0.03]);
xlim([-pi pi])


if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-DecoderErrorM40g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-DecoderErrorM40g0_5rnd1N2k.pdf')
    end
end

figure;plot(angle(result.with_learning{1}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])
if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM6g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM6g0_5rnd1N2k.pdf')
    end
end
figure;plot(angle(result.with_learning{2}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])
if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM12g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM12g0_5rnd1N2k.pdf')
    end
end
figure;plot(angle(result.with_learning{3}.z_cl(1:end,:))','.r');
box off;ylim([-pi pi]);xlim([0 1000])

if saveFlag
    if net.g==1.5
        saveas(gcf,'FigSampling-PsiVsTimeM40g1_5rnd1N2k.pdf')
    elseif net.g==0.5
        saveas(gcf,'FigSampling-PsiVsTimeM40g0_5rnd1N2k.pdf')
    end
end

