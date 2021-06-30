%%  Ring figure
% first run nn_different_rotations_taueff_Delta with:
% DegIC_vec=[270,300,0, 180,120]; FinalFlag=1; DegRot_vec=[90 90 90 90 90]; h_vec=0; g_vec=1;
% load('FigureResultsN4000a2Deg90h0g1eps0p01seed1.mat')
st=50;
jump=1;
ColorsMat=hsv;
figure
for ddiicc=1
    % figure
    hh=plot(zdeltaRot{ddiicc}(st:jump:end),'o','MarkerSize',5,'Color',ColorsMat(10*(ddiicc),:));
    hold on;
    plot(sim.f_ol(1,:),sim.f_ol(2,:),'-k');
end
box off; axis square; xlim([-1.5 1.5])

figure

for ddiicc=1
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
g_vec=[0.01 1 1.5]
for g=g_vec
%     RotDeg
%     Rot90DegVsg_ResultsN1000a2Deg90h0g0p01eps0p01seed1
    load(['Rot90DegVsg_ResultsN1000a2Deg90h0','g',strrep(num2str(g),'.','p'),'eps0p01seed1.mat'])
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
    axis square; box off;  xlim([0 1000]);ylim([0 pi/2])
end

