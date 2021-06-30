% Do the closed loop
phi=@(x)erf(x./sqrt(2));
PlotFlag=0;
PlotEigFlag=1;
% K_vec=[2:2:30];
K_vec=[1,2,4,6, 20];
% for ii=[1 3 4]
for ii=[2]
    switch ii
        case 1
            N=1000;a=4;A=1.2; g=0.2; seed=1; UniformFlag=0;M=40;h=0.5; % classical ring
        case 2
            N=4000;a=4;A=0.9; g=0.2; seed=1; UniformFlag=0;M=40;h=-0.4; % Heterogenous ring
    end
end

% [result]=nn_plot_tuning_curves(a,h,A,g,N,M,seed,UniformFlag,PlotFlag,PlotEigFlag);
[result]=nn_plot_compare_to_semiemp(a,h,A,g,N,M,seed,UniformFlag,PlotFlag,PlotEigFlag,K_vec);

%%

% r=(phi(result.with_learning.x_cl));
% r=(phi(result.with_learning.x_cl));
r=(phi(result.sim.x));

Crr=r*r';  [v,lambda]=eig(Crr);
v=1/length(Crr)*1/N*v; %doublecheck normalization

[lambda,ii]=sort(diag(lambda),'descend');
lambda=abs(real(lambda))+1e-20;
PR=sum(lambda)^2/sum(lambda.^2)

v=v(:,ii); v0=diag(v(:,1));
figure(10)
semilogy(lambda,'o');xlim ([0 M])
% ylim([1e-18 1e0]);
box off;axis square;hold on;

% CrrT=r'*r./N;

% [eta,lambda]=eig(CrrT);
% [lambda,ii]=sort(diag(lambda),'descend');
% plot(lambda)

%%
maxPlotEvR=8;maxPlotEvPsi=8;
    figure
    kk=4;
    %     K_vec=[2:2:24];
    %     Kmax=K_vec(kk);
    
    plot(real(result.ee_cl1(3:end)),imag(result.ee_cl1(3:end)),'.y')
    
    hold on
    plot(real(result.ee_cl1(1:maxPlotEvR)),imag(result.ee_cl1(1:maxPlotEvR)),'.g')
    
    plot(real(result.ee_cl2(1:maxPlotEvPsi)),imag(result.ee_cl2(1:maxPlotEvPsi)),'.r')
    
    plot( real(result.y{kk,1}),imag(result.y{kk,1}),'og','MarkerSize',10)
    hold on
    plot( real(result.y{kk,2}),imag(result.y{kk,2}),'or','MarkerSize',10)
    plot( real(result.y_all{kk}),imag(result.y_all{kk}),'*m','MarkerSize',10)
    
    r=1+result.bulkMF;x=-1;y=0; th = 0:pi/60:2*pi;xunit = r * cos(th) + x; yunit = r * sin(th) + y;
    plot(xunit, yunit,'--');
    
    box off; axis square
%     set(gcf,'Position',[   977   748   143   200])
    
    
