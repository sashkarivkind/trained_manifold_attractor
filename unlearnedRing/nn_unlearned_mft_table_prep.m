clear;
tmax=1000;

N=4000;
net.phi=@(x)erf(x./sqrt(2)); 
net.phip =@(x) exp(-x.^2/2)/sqrt(pi/2);

% g_vec=[0.01, 0.1, 0.3, 0.5];
% A_vec=[0.5,1,1.2,1.5,2.];

g_vec=[0.01, 0.1:0.1:2];
A_vec=[0.5,1,1.2,1.5,2.];
%% cover up for forgetting to save J2
uu = load('A_J2_table');
J2_vec=[];
for aa = 1:length(A_vec)
 ii=find(uu.A_vec==A_vec(aa));
 if length(ii)~=1
     error
 end
 J2_vec(ii)=uu.J2_vec(ii);
end
%%

seed_vec=1:5;
rmse_rec=[];
A_rec=[];
g_rec=[];
seed_rec=[];
cnt=0;
for aa=1:length(A_vec)
    A=A_vec(aa);
    for gg=1:length(g_vec)
        g=g_vec(gg);
        uu=struct;
[uu.mu,uu.sig,uu.A,uu.d,uu.d_approx]=solveFixedPointRingUnlearned(net.phi,0,g,trace(J2_vec(aa)),0.,1.,0,0,net.phip);
    mft_result_table{aa,gg}=uu;
    uu
    end
end
