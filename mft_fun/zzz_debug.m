sig=1e-0;
A=1.2;
mu=0;
D=@(z)exp(-((z.^2)/2))./(sqrt(2*pi));
dtheta=1./(2*pi);
phi=net.phi;

m2 =@(theta) arrayfun(@(theta) (sin(theta) ...
    .*integral(@(z) phi(mu+sig.*z+A*cos(theta)).*D(z),-Inf,Inf)).^2.*dtheta,theta);
m2int= integral( m2, 0,2*pi);
v2=@(theta,z) (sin(theta).*phi(mu+sig.*z+A*cos(theta))).^2.*D(z).*dtheta;
v2int=integral2(v2,0,2*pi,-Inf,Inf)
d2=v2int-m2int;
d=sqrt(d2);

% %%
% disp('----------')
% theta=1.3
% (sin(theta) ...
%     .*integral(@(z) phi(mu+sig.*z+A*cos(theta)).*D(z),-Inf,Inf)).^2.*dtheta
% 
% (sin(theta).*phi(mu+sig.*z+A*cos(theta))).^2.*dtheta
% 
% 
% %%
% disp('----------')
% (integral(@(z) phi(sig.*z+A*cos(theta)).*D(z),-Inf,Inf))
% 
% (phi(sig.*z+A*cos(theta)))