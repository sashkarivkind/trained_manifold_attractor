function [mu_out, sig_out,A_out,d_out,d_out_approx]=solveFixedPointRingUnlearned(phi,I,J_0,J_2,mu_0,sig_0,K,Delta_vec,phip)
% Ver1: 27/07/18
% Function solves the fixed point statistics of a ring atracttor of
% given J_2.
% with erf: r=phi(x)=erf(x/sqrt(2))=int_0^xdye^(-y^2/2)1/sqrt(2pi)?
%
%
% Inputs:
%       phi             Transfer function
%       I               External sqrt(K) input
%       J_0             Synaptic strength stochastic
%       J_2             Synaptic strength structured
%       mu_0,sig_0      initial conditions for mu and sigma
%       K               num of inputs
%       Delta_vec       angles to be learned
% Outputs
%       mu_out          mean of the field
%       sig_out         sd of the field
%       beta_out        correlations of the field. beta_out(theta=0)=sig_out^2
%       A               Amplitude of the bump

D=@(z)exp(-((z.^2)/2))./(sqrt(2*pi));
dtheta=1./(2*pi);
r=I/J_0;

options=optimset('Display','on');% Tae off if ou want to see display!

%Solve mu an sigma with arbitrary TF
A_0=J_2;
Abias=1e-5;
a0=[sig_0,A_0];
y=fsolve(@F,a0,options);
mu_out=0;sig_out=y(1); A_out=y(2)+Abias;
[d_out,d_out_approx]=RmsDeltaParallel(sig_out,A_out);
% %Solve beta (corrleatoons in field) with Relu
% ind=0;
% a0=sig_out.^2;
%
% y2=zeros(size(Delta_vec));
% beta_add=zeros(size(Delta_vec));
%
% for Delta=Delta_vec
%     ind=ind+1;
% %     tic;
%     y2(ind)=fsolve(@G_1Erf,a0,options);
% %     toc;
%     a0=y2(ind)+1e-3;
%     beta_add(ind)=integral(@(theta)cos(theta).*cos(theta+Delta).*dtheta,0,2*pi);
% end
% beta_out=y2;
%
% % Add the variance due to the FF input
% sig_out=sqrt(sig_out.^2+A^2/2);
% beta_out=beta_out+A.^2.*beta_add;


    function y=F(a)
        mu=0;
        sig=a(1);
        A=Abias+a(2);
        %         g=@(theta,z) phi(mu+sig.*z+A*cos(theta)).*D(z).*dtheta;
        g2=@(theta,z) phi(mu+sig.*z+A*cos(theta)).^2.*D(z).*dtheta;
        m2=@(theta,z) cos(theta).*phi(mu+sig.*z+A*cos(theta)).*D(z).*dtheta;
        %         y(1)=integral2(g,0,2*pi,-Inf,Inf)-r+mu/sqrt(K);
        %         y(1)=integral2(g,0,2*pi,-Inf,Inf)-r;
        y(1)=sqrt(J_0^2*integral2(g2,0,2*pi,-Inf,Inf));
        y(2)=J_2*integral2(m2,0,2*pi,-Inf,Inf)-Abias;
        y=y-a;
    end

    function [d,d_approx]=RmsDeltaParallel(sig,A)
        mu=0;
        m2 =@(theta) arrayfun(@(theta) (J_2*sin(theta) ...
            .*integral(@(z) phi(mu+sig.*z+A*cos(theta)).*D(z),-Inf,Inf)).^2.*dtheta,theta);
        m2int= integral( m2, 0,2*pi);
        v2=@(theta,z) (J_2*sin(theta).*phi(mu+sig.*z+A*cos(theta))).^2.*D(z).*dtheta;
        v2int=integral2(v2,0,2*pi,-Inf,Inf);
        d2=v2int-m2int;
        d=sqrt(d2);
        v2_approx=@(theta) (J_2*sig*sin(theta).*phip(mu+A*cos(theta))).^2.*dtheta;
        d2_approx=integral(v2_approx,0,2*pi);
        d_approx=sqrt(d2_approx);
    end

    function y=G_1(a)
        beta=a;
        a1=@(theta,z3) mu_out+sign(beta).*sqrt(abs(beta)).*z3+A*cos(theta);
        a2=@(theta,z3)mu_out+sqrt(abs(beta)).*z3+A*cos(theta+Delta);
        b=sqrt(sig_out.^2-abs(beta));
        if b<1e-9
            b=1e-8;
        end
        
        F1a= @(theta,z3) a1(theta,z3)./2.*(      1+erf( a1(theta,z3)./sqrt(2.*b)  )      ) +b./sqrt(2*pi).*exp(-a1(theta,z3).^2./(2.*b.^2));
        F1b= @(theta,z3) a2(theta,z3)./2.*(      1+erf( a2(theta,z3)./sqrt(2.*b)  )      ) +b./sqrt(2*pi).*exp(-a2(theta,z3).^2./(2.*b.^2));
        
        g4=@(theta,z3) F1a(theta,z3).*F1b(theta,z3).*D(z3).*dtheta;
        gall=integral2(g4,0,2*pi,-Inf,Inf);
        y=J_0^2.*gall-beta;
    end

    function y=G_1Erf(a)
        beta=a;
        a1=@(theta,z3) mu_out+sign(beta).*sqrt(abs(beta)).*z3+A*cos(theta);
        a2=@(theta,z3)mu_out+sqrt(abs(beta)).*z3+A*cos(theta+Delta);
        b=sqrt(sig_out.^2-abs(beta));
        if b<1e-9
            b=1e-8;
        end
        
        F1a= @(theta,z3) erf(a1(theta,z3)./sqrt(2*(1+b.^2)));
        F1b= @(theta,z3) erf(a2(theta,z3)./sqrt(2*(1+b.^2)));
        
        g4=@(theta,z3) F1a(theta,z3).*F1b(theta,z3).*D(z3).*dtheta;
        gall=integral2(g4,0,2*pi,-Inf,Inf);
        y=J_0^2.*gall-beta;
    end


%     function y=G(a)
%         beta=a;
%         g4=@(theta,z1,z2,z3) phi(mu_out+sqrt(sig_out.^2-abs(beta)).*z1+sign(beta).*sqrt(abs(beta)).*z3+A*cos(theta)).*phi(mu_out+sqrt(sig_out.^2-abs(beta)).*z2+sqrt(sig_out.^2-abs(beta)).*z3+A*cos(theta+Delta)).*D(z1).*D(z2).*D(z3).*dtheta;
%         gall=integral(@(theta)integral3(@(z1,z2,z3)g4(theta,z1,z2,z3),-Inf,Inf,-Inf,Inf,-Inf,Inf),0,2*pi,'ArrayValued',true);
%         %     g4=@(theta,z1,z2,z3) phi(mu_outs+sqrt(sig_out.^2-abs(beta)).*z1+sign(beta).*sqrt(abs(beta)).*z3+A*cos(theta)).*phi(mu_out+sqrt(sig_out.^2-abs(beta)).*z2+sqrt(sig_out.^2-abs(beta)).*z3+A*cos(theta+Delta)).*D(z1).*D(z2).*D(z3).*dtheta;
%         y=J_0^2.*gall-beta;
%     end



end





