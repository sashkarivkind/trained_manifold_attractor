figure;
for uu=1:10
plot(phase([1,1i]*result.with_learning{3,2,uu}.z_ol)-sim.psi);
hold on;

end
%%
mm=5;
xx=[];
yy=[];
for uu=1:10
v_ol=phase([1,1i]*result.with_learning{3,mm,uu}.z_ol)-sim.psi;
v_cl=diff(mod(angle(result.with_learning{3,mm,uu}.z_cl.'),2*pi))/dt;
xx(:,end+1)=v_ol;
yy(:,end+1)=v_cl(end,:);
end
figure;
plot(xx(:,:),yy(:,:),'x');


%% just plot drift velocities
mm=1;
smpl_pnt_num=3;
xx=[];
yy=[];
figure;
for gg=2
for uu=1:10
v_ol=phase([1,1i]*result.with_learning{gg,mm,uu}.z_ol)-sim.psi;
v_cl=diff(mod(angle(result.with_learning{gg,mm,uu}.z_cl.'),2*pi))/dt;

v_cl=v_cl(:,2:end);% xx(:,end+1)=v_ol;
plot(v_cl');
hold on;
yy(gg,uu)=min((v_cl(smpl_pnt_num,:)));
end
end
figure;
plot(g_vec,mean(yy(:,:),2),'x-');