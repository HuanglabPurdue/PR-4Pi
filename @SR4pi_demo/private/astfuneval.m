function [SSE,Sx,Sy] = astfuneval(params,z,inSx,inSy,fitType)
Ax = params(1);
Bx = params(2);
Ay = params(3);
By = params(4);
gamma = params(5);
d = params(6);
w = params(7);
Sx = w*sqrt(1+((z-gamma)/d).^2+Ax.*((z-gamma)/d).^3+Bx.*((z-gamma)/d).^4);
Sy = w*sqrt(1+((z-gamma)/d).^2+Ay.*((z-gamma)/d).^3+By.*((z-gamma)/d).^4);
SSE = 0;
if nargin>2
    SSE1=sum((inSx-Sx).^2);
    SSE2=sum((inSy-Sy).^2);
    switch fitType
        case 1
            SSE=SSE1+SSE2;
        case 2
            SSE=SSE1;
        case 3
            SSE=SSE2;
    end
end
end