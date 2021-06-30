%%  Ring figure
% first run nn_different_rotations_taueff_Delta with:
% DegIC_vec=[270,300,0, 180,120]; FinalFlag=1; DegRot_vec=[90 90 90 90 90]; h_vec=0; g_vec=1;
load('FigureResultsN4000a2Deg90h0g1eps0p01seed1.mat')
st=50;
jump=1;
ColorsMat=hsv;
figure
for ddiicc=1:length(DegIC_vec)
    % figure
    hh=plot(zdeltaRot{ddiicc}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(ddiicc),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
end
box off; axis square; xlim([-1.5 1.5])

figure

for ddiicc=1:size(zdeltaRot,2)
    psi_dot=diff(angle(zdeltaRot{ddiicc}))/dt;
    
    hhh=plot(angle(zdeltaRot{ddiicc}(st-1+2:jump:end)),psi_dot(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(ddiicc),:))
    hold on
end
psi_on_manifold = (atan2(sim.f_ol(2,:),sim.f_ol(1,:)));
Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
plot_with_sorted_x(psi_on_manifold,Delta./tau_eff2,'k-')
plot_with_sorted_x(psi_on_manifold,Delta,'k:')
box off; axis square;

xlim([-pi pi]);ylim([-1 1]*0.01);set(gca,'Xtick',[-pi:pi/4:pi])
%% Ellipse figure
% first run nn_different_rotations_taueff_Delta with:
%h_vec=0.2; DegRot=180; DegIC_vec=[0,90]; FinalFlag=0; DegRot_vec=[182,180];

ind=0;
ColorsMat=hsv;
for RotDeg=[182 180 ]
    RotDeg
    load(['FigureResultsN4000a2Deg', num2str(RotDeg) ,'h0p2g1eps0p01seed1.mat'])
    ind=ind+1
    st=50;
    jump=1;
    
    
    figure
    psi_dot=diff(angle(zdeltaRot{1}))/dt;
    
    hhh=plot(angle(zdeltaRot{1}(st-1+2:jump:end)),psi_dot(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(ind),:));
    hold on
    psi_on_manifold = (atan2(sim.f_ol(2,:),sim.f_ol(1,:)));
    Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
    plot_with_sorted_x(psi_on_manifold,Delta./tau_eff2,'k-')
    plot_with_sorted_x(psi_on_manifold,Delta,'k:')
    box off; axis square;
    
    xlim([-pi pi]);ylim([-1 1]*0.015);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    
    
    figure(9);
    hh=plot(zdeltaRot{1}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(ind),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
    
    box off; axis square; xlim([-1.5 1.5])
    
    figure(10);
    plot(abs(unwrap(((angle(zdeltaRot{1})-angle(zdeltaRot{1}(1)))))),'Color',ColorsMat(10*(ind),:));
    hold on;
end

%% yppo figure
% first run nn_different_rotations_taueff_Delta with:
%h_vec=0.2; DegRot=180; DegIC_vec=[0,90]; FinalFlag=0; DegRot_vec=[182,180];

jjj=0;
ColorsMat=hsv;
for RotDeg=[135 136 ]
    RotDeg
    load(['FigureResultsN4000a4Deg', num2str(RotDeg) ,'h-0p15g1eps0p01seed1.mat'])
    
    jjj=jjj+1
    st=50;
    jump=1;
    
    
    figure
    psi_dot=diff(angle(zdeltaRot{1}))/dt;
    
    hhh=plot(angle(zdeltaRot{1}(st-1+2:jump:end)),psi_dot(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on
    psi_on_manifold = (atan2(sim.f_ol(2,:),sim.f_ol(1,:)));
    Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
    plot_with_sorted_x(psi_on_manifold,Delta./tau_eff2,'k-')
    plot_with_sorted_x(psi_on_manifold,Delta,'k:')
    box off; axis square;
    
    xlim([-pi pi]);ylim([-1 1]*0.01);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    
    
    figure(9);
    hh=plot(zdeltaRot{1}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
    
    box off; axis square; xlim([-1.5 1.5])
    
    figure(10);
    plot(abs(unwrap(((angle(zdeltaRot{1})-angle(zdeltaRot{1}(1)))))),'Color',ColorsMat(10*(jjj),:));
    hold on;
end

%% Not so Bad Hyppo figure
% first run nn_different_rotations_taueff_Delta with:
%h_vec=0.2; DegRot=180; DegIC_vec=[0,90]; FinalFlag=0; DegRot_vec=[182,180];

jjj=0;
ColorsMat=hsv;
for RotDeg=[90 91 ]
    RotDeg
    load(['FigureResultsN4000a4Deg', num2str(RotDeg) ,'h-0p2g1eps0p01seed1.mat'])
    
    jjj=jjj+1
    st=50;
    jump=1;
    
    
    figure
    psi_dot=diff(angle(zdeltaRot{1}))/dt;
    
    hhh=plot(angle(zdeltaRot{1}(st-1+2:jump:end)),psi_dot(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on
    psi_on_manifold = (atan2(sim.f_ol(2,:),sim.f_ol(1,:)));
    Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
    plot_with_sorted_x(psi_on_manifold,Delta./tau_eff2,'k-')
    plot_with_sorted_x(psi_on_manifold,Delta,'k:')
    box off; axis square;
    
    xlim([-pi pi]);ylim([-1 1]*0.01);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    
    
    figure(9);
    hh=plot(zdeltaRot{1}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
    
    box off; axis square; xlim([-1.5 1.5])
    
    figure(10);
    plot(abs(unwrap(((angle(zdeltaRot{1})-angle(zdeltaRot{1}(1)))))),'Color',ColorsMat(10*(jjj),:));
    hold on; box off; axis square;
    set(gca,'ytick',[0 pi/4, pi/2])
end

%% Very Bad Hyppo figure
% first run nn_different_rotations_taueff_Delta with:
%h_vec=0.2; DegRot=180; DegIC_vec=[0,90]; FinalFlag=0; DegRot_vec=[182,180];

jjj=0;
ColorsMat=hsv;
for RotDeg=[90 91 160 ]
    RotDeg
    load(['FigureResultsN4000a4Deg', num2str(RotDeg) ,'h0p9g1eps0p01seed1.mat'])
    
    jjj=jjj+1
    st=50;
    jump=1;
    
    
    figure
    psi_dot=diff(angle(zdeltaRot{1}))/dt;
    
    hhh=plot(angle(zdeltaRot{1}(st-1+2:jump:end)),psi_dot(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on
    psi_on_manifold = (atan2(sim.f_ol(2,:),sim.f_ol(1,:)));
    Delta=angle(z_epsi_on_manifold(1,:)+1i*z_epsi_on_manifold(2,:))-psi_on_manifold;
    plot_with_sorted_x(psi_on_manifold,Delta./tau_eff2,'k-')
    plot_with_sorted_x(psi_on_manifold,Delta,'k:')
    box off; axis square;
    
    xlim([-pi pi]);ylim([-1 1]*0.01);set(gca,'Xtick',[-pi:pi/4:pi])
    
    
    
    
    figure(9);
    hh=plot(zdeltaRot{1}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(jjj),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
    
    box off; axis square; xlim([-1.5 1.5])
    
    figure(10);
    plot(abs(unwrap(((angle(zdeltaRot{1})-angle(zdeltaRot{1}(1)))))),'Color',ColorsMat(10*(jjj),:));
    hold on; box off; axis square;
    set(gca,'ytick',[0 pi/4, pi/2])
end
%% Plot EV
RotDeg=[90]
load(['FigureResultsN4000a4Deg', num2str(RotDeg) ,'h0p9g1eps0p01seed1.mat'])
% load(['FigureResultsN4000a4Deg', num2str(RotDeg) ,'h-0p2g1eps0p01seed1.mat'])
dd=1;
ind=0;
for pptt=[1  20]
    ind=ind+1;
% pptt=1;
Sfull=(net.W+net.wfb*net.wout')*diag(net.phip(sim.x(:,pptt)))-eye(net.N);
S1=(net.W+net.wfb(:,1)*net.wout(:,1)')*diag(net.phip(sim.x(:,pptt)))-eye(net.N);
S2=(net.W+net.wfb(:,2)*net.wout(:,2)')*diag(net.phip(sim.x(:,pptt)))-eye(net.N);
eSfull=eig(Sfull);
ev_full_vec{ind}=sort(eSfull);
eS1=eig(S1);
ev_full_d1{ind}=sort(eS1);
eS2=eig(S2);
ev_full_d2{ind}=sort(eS2);

subplot(2,2,2*ind-1)
plot(real(eS1),imag(eS1),'.g')
hold on
plot(real(eS2),imag(eS2),'.r')
title(sprintf('Rho: Point=%d, MaxEv inR=%.4f',pptt,max(real(eS1)) ))

xlim([-1.5 0.1])

subplot(2,2,2*ind-1+1)
plot(real(eSfull),imag(eSfull),'.')
title(sprintf(' Full: Point=%d, MaxEv inR=%.4f',pptt,max(real(eSfull)) ))
xlim([-1.5 0.1])

end
% plot(real(eS),imag(eS),'.')