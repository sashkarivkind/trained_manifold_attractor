g_vec=[0.01 0.1:0.1:1.8];
% o = get_betas_etc( g_vec,'pp_p1' );
o = get_betas_etc( g_vec,'betapp128' );

for gg=1:length(g_vec)
% beta(gg)=(1-o{gg}.v2(1,1));
beta(gg)=(o{gg}(1,1));
end
plot(g_vec,beta); hold on;
axis square
xlim([0 1.8]);ylim([0 0.7]);box off;
xlabel('g'); ylabel('beta_1(1,1)');title('MF with hmax=5')
%%
g_vec=[0.01 0.1:0.1:1.8];
% o = get_betas_etc( g_vec,'pp_p1' );
o = get_betas_etc( g_vec,'betapp' );

for gg=1:length(g_vec)
% beta(gg)=(1-o{gg}.v2(1,1));
beta(gg)=(o{gg}{2}(1,1));
end
    plot(g_vec,beta)
hold on;
%%
o = get_betas_etc( g_vec,'pp_p1' );

for gg=1:length(g_vec)
beta(gg)=(1-o{gg}.v2(1,1));
end
    plot(g_vec,beta)
