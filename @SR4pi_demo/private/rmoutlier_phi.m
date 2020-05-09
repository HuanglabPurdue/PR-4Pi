function [mask,zT,phi0] = rmoutlier_phi(phi1,z1,phic,plotflag)

f = @(x,z,phi) sum((wrapToPi(z.*x(1)+x(2))-phi).^2);
phi = phi1;
z = z1;
[est,sse] = fminsearch(@(x) f(x,z,phi),[2*pi/0.28,wrapToPi(phic)]);

y = wrapToPi(z.*est(1)+est(2));
zT = 2*pi/est(1)*1e3;
phi0 = est(2);
phithresh = 40*pi/zT;

mask = abs(wrapToPi(phi-y))<phithresh;

if plotflag == 1
    
    figure;
    plot(z,phi,'.',z(mask),phi(mask),'.')
    figure;
    plot(z,wrapToPi(phi-y),'.',z(mask),wrapToPi(phi(mask)-y(mask)),'.')
    
end
