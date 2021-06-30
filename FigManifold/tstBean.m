% Do the closed loop
phi=@(x)erf(x./sqrt(2));
PlotFlag=1;
PlotEigFlag=1;
% for ii=[1 3 4]

N=1000;a=0.5;A=1.2; g=1; seed=1; UniformFlag=0;M=40;h=0; % classical ring



[result]=nn_plot_tuning_curves_Bean(a,h,A,g,N,M,seed,UniformFlag,PlotFlag,PlotEigFlag);
% Plot the tuning curves
%%
id=[1:250:N];
xaxis_deg=(0:(size(result.with_learning.z_ol,2)-1))*2*pi/size(result.with_learning.z_ol,2);
figure
plot(xaxis_deg-pi,(phi(result.with_learning.x_cl(id,:)))');
xlim([-pi pi]); ylim([-1 1]); box off ;axis square;set(gca,'Xtick',[-pi:pi/4:pi])
title(['a= ',num2str(a), ', h= ',num2str(h), ', g= ',num2str(g) ])

%%

% Plot the correaltion function
% inds=randperm(1000);
r=(phi(result.with_learning.x_cl));
% r=r(inds,:);

 Crr=r*r';  [v,lambda]=eig(Crr);
                v=1/length(Crr)*1/N*v; %doublecheck normalization
                
                [lambda,ii]=sort(diag(lambda),'descend');
                lambda=abs(real(lambda))+1e-20;
                v=v(:,ii); v0=diag(v(:,1));
figure(10)
semilogy(lambda,'o');xlim ([0 M])
% ylim([1e-18 1e0]);
box off;axis square;hold on;


% plot(cumsum(flipud(diag(D)))./sum(flipud(diag(D))))
%
%%
% % % figure
% % % subplot 311
% % % plot(lambda(1).*v(:,1)'*r,lambda(2).*v(:,2)'*r,'o')
% % % ylim([min(lambda(1).*v(:,1)'*r) max(lambda(1).*v(:,1)'*r)])
% % % xlim([min(lambda(1).*v(:,1)'*r) max(lambda(1).*v(:,1)'*r)])
% % % axis square; box off
% % % subplot 312
% % % % plot3(lambda(1).*v(:,1)'*r,lambda(2).*v(:,2)'*r,lambda(7).*v(:,7)'*r,'o')
% % % plot(lambda(3).*v(:,3)'*r,lambda(4).*v(:,4)'*r,'o')
% % % ylim([min(lambda(3).*v(:,3)'*r) max(lambda(3).*v(:,3)'*r)])
% % % xlim([min(lambda(3).*v(:,3)'*r) max(lambda(3).*v(:,3)'*r)])
% % % axis square; box off
% % % 
% % % subplot 313
% % % plot(lambda(5).*v(:,5)'*r,lambda(6).*v(:,6)'*r,'o')
% % % ylim([min(lambda(5).*v(:,5)'*r) max(lambda(5).*v(:,5)'*r)])
% % % xlim([min(lambda(5).*v(:,5)'*r) max(lambda(5).*v(:,5)'*r)])
% % % axis square; box off
% % % 
% % % figure
% % % % subplot 311
% % % plot3(lambda(1).*v(:,1)'*r,lambda(2).*v(:,2)'*r,lambda(3).*v(:,3)'*r,'o')
% % % ylim([min(lambda(1).*v(:,1)'*r) max(lambda(1).*v(:,1)'*r)])
% % % xlim([min(lambda(1).*v(:,1)'*r) max(lambda(1).*v(:,1)'*r)])
% % % 
% % % zlim([min(lambda(1).*v(:,1)'*r) max(lambda(1).*v(:,1)'*r)])
% % % axis square; box off

figure
subplot 311
plot(v(:,1)'*r,v(:,2)'*r,'o')
ylim([min(v(:,1)'*r) max(v(:,1)'*r)])
xlim([min(v(:,1)'*r) max(v(:,1)'*r)])
axis square; box off
subplot 312
% plot3(lambda(1).*v(:,1)'*r,lambda(2).*v(:,2)'*r,lambda(7).*v(:,7)'*r,'o')
plot(v(:,3)'*r,v(:,4)'*r,'o')
ylim([min(v(:,3)'*r) max(v(:,3)'*r)])
xlim([min(v(:,3)'*r) max(v(:,3)'*r)])
axis square; box off

subplot 313
plot(v(:,5)'*r,v(:,6)'*r,'o')
ylim([min(v(:,5)'*r) max(v(:,5)'*r)])
xlim([min(v(:,5)'*r) max(v(:,5)'*r)])
axis square; box off

figure
% subplot 311
colorVec = hsv(size(r,2));
for ll=1:1:size(r,2)
ll
plot3(v(:,1)'*r(:,ll),v(:,2)'*r(:,ll),v(:,3)'*r(:,ll),'o','Color',colorVec(ll,:));
hold on
end
% % plot3(v(:,1)'*r,v(:,2)'*r,v(:,3)'*r,'o')
% % hold on
plot3(v(:,1)'*r,v(:,2)'*r,min(v(:,1)'*r).*ones(size(v(:,3)'*r))./2,'o')
grid on

ylim([min(v(:,1)'*r) max(v(:,1)'*r)])
xlim([min(v(:,1)'*r) max(v(:,1)'*r)])

zlim([min(v(:,1)'*r) max(v(:,1)'*r)]./2)
axis square; box off

%%
% maxPlotEvR=4;maxPlotEvPsi=2;

maxPlotEvR=8;maxPlotEvPsi=8;
if PlotEigFlag
    figure
    kk=6;
%     K_vec=[2:2:24];
%     Kmax=K_vec(kk);
    
    plot(real(result.ee_cl1(3:end)),imag(result.ee_cl1(3:end)),'.y')
    
    hold on
    plot(real(result.ee_cl1(1:maxPlotEvR)),imag(result.ee_cl1(1:maxPlotEvR)),'.g')
    
    plot(real(result.ee_cl2(1:maxPlotEvPsi)),imag(result.ee_cl2(1:maxPlotEvPsi)),'.r')
    
    plot( real(result.y{kk,1}),imag(result.y{kk,1}),'og','MarkerSize',10)
    hold on
    plot( real(result.y{kk,2}),imag(result.y{kk,2}),'or','MarkerSize',10)
    
        r=1+result.bulkMF;x=-1;y=0; th = 0:pi/60:2*pi;xunit = r * cos(th) + x; yunit = r * sin(th) + y;
        plot(xunit, yunit,'--');
        
%     xlim ([-2  0.01]); ylim([-0.3 0.3]); box off; axis square
%     xlim ([-1.2  0.04]); ylim([-0.3 0.3]); 
% xlim ([-1.5  0.01]); ylim([-0.15 0.15]); 
% xlim ([-1.5  1]); ylim([-0.15 0.15])
box off; axis square
set(gcf,'Position',[   977   748   143   200])
    
end
% %%
% r=(phi(result.with_learning.x_cl));
%  Crr=r'*r;  [eta,lambda]=eig(Crr);
%                 v=1/length(Crr)*1/N*v; %doublecheck normalization
%                 
%                 [lambda,ii]=sort(diag(lambda),'descend');
%                 lambda=abs(real(lambda))+1e-20;
%                 eta=eta(:,ii); 