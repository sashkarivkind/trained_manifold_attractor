function crrfu=crr_fourier_fun(sigma,g,A,phi,K,nmc)

if nargin < 6
    nmc=1e7;
end
%%

uu1=rand(nmc,2);
zz1=randn(nmc,K);
zz2=randn(nmc,K);

% %%
% u=0.1:0.01:2;
% v=zeros(size(u));
% for uu=1:length(u)
%     v(uu)=c_iter(u(uu));
% end
% figure;
% plot(u,v);
%%
% options = optimoptions('fsolve','Display','iter');
c_ini=rand(K,1);
crrfu=fsolve(@(x) c_iter(x) - x ,c_ini);


    function CC_=c_iter(CC)
        tp=uu1*2*pi; % distribution of angles theta and psi
        mm=real(sqrt(CC));
        ai=zz1*diag(mm); %cos coefficient
        bi=zz2*diag(mm); %sin coefficient
        harmonics=tp(:,2)*[(1:length(CC))*2-1]; %taking odd coefficients only
        x1_approx=sum((ai.*cos(harmonics))...
            +(bi.*sin(harmonics)),2);
        
        CC_=zeros(length(CC),1);
        for ii=1:length(CC)
            CC_(ii)=g^2*mean(...
                cos((2*ii-1)*tp(:,2)).*...
                phi(A*cos(tp(:,1)-tp(:,2))+x1_approx)...%phi integrated over psi
                .*phi(A*cos(tp(:,1))+sum(ai,2))...%at origin we have a sum of all the cosine coeffitients
                );
            
        end
        [CC,CC_]
    end
end
%% sdsds
%     0.9597
%     0.0428
%     0.0072
%     0.0024
%    -0.0001

%    0.9618
%     0.0444
%     0.0044
%     0.0000
%    -0.0000

%0.4272    0.0183    0.0025    0.0005    0.0002
%0.8544    0.0366    0.0051    0.0009    0.0004
% 
% crrfu =
% 
%     0.4258
%     0.0088
%     0.0004
%    -0.0003
%    -0.0004
% 
% crrfu =
% 
%     0.4259
%     0.0069
%     0.0005
%    -0.0000
%    -0.0003
