function w=lms_weights(r,f,regfac)
if nargin<3
    regfac=0;
end
w=(r/(r'*r+regfac))*f';