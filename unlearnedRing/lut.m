function ii=lut(vec,val,en_mult_val,epsilon)
if nargin<3
    en_mult_val=0;
end
if nargin<4
    epsilon=1e-10;
end
ii=find(abs(vec-val)<epsilon);
if ~en_mult_val&&(length(ii)~=1)
    error
end