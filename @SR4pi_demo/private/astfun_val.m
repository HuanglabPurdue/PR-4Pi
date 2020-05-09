function [FittedCurve] = astfun_val(zdata,params)
w = params(1);
c = params(2);
d = params(3);
A = params(4);
B = params(5);
FittedCurve = w.*sqrt(1+((zdata-c)./d).^2+A.*((zdata-c)./d).^3+B.*((zdata-c)./d).^4);
end