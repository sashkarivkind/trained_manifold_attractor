function o=my_prep_errb(x,y)
o.xx=sort(unique(x));
for jj=1:length(o.xx)
    this_x=o.xx(jj);
    yy=y(x==this_x);
    o.yy_m(jj)=mean(yy);
    o.yy_std(jj)=std(yy);
    o.yy_n(jj)=length(yy);
end
