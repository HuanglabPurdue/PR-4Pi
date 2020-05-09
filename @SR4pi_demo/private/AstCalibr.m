function [SSE,Sx,Sy]=AstCalibr(x,inSx,inSy,z,fitType)
Ax=x(1);               %Aberration Terms
Bx=x(2);
Ay=x(3);
By=x(4);
gamma=x(5) ;          %separation of x/y focal planes
d=x(6);
PSFsigmax=x(7);
PSFsigmay=x(8);
z0 = x(9);
Sx=PSFsigmax*sqrt(1+((z-z0-gamma)/d).^2+Ax*((z-z0-gamma)/d).^3+Bx*((z-z0-gamma)/d).^4);
Sy=PSFsigmay*sqrt(1+((z-z0+gamma)/d).^2+Ay*((z-z0+gamma)/d).^3+By*((z-z0+gamma)/d).^4);

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
