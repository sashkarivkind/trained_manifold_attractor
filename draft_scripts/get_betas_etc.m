function o = get_betas_etc( g_vec,field )
%GET_BETAS_ETC Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    field='betapp';
end
beta_info=load('beta_graph_long.mat');
for gg=1:length(g_vec)
    ii=find(abs(beta_info.J_vec-g_vec(gg))<1e-10);
    if length(ii)~=1
        error('pole info not found!');
    end
    o{gg}=beta_info.debuu{ii}.(field);
end

end

