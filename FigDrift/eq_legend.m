function o=eq_legend(vec,pref)
o={};
for vv=1:length(vec)
    o{end+1}=[pref,num2str(vec(vv))];
end