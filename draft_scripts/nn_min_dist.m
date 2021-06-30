function o=nn_min_dist(v1,v2)
%for v1 returns minimal distances from any element of v2

V1=v1*ones(size(v1))';
V2=v2*ones(size(v2))';
o=min(abs(V1-V2.'),[],2);