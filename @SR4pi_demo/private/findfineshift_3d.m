function [xshift, yshift, zshift,xout, yout,zout] = findfineshift_3d(x,y,z,cormask,sigma,stepsize)
%findfineshift Auto Correct Drift in SR_demo object
%   Cross correlation between image stacks.
% INPUT
%   cormask: N x 1 integer array grouping observations for corellation (eg filenumber)
%   sigma: smoothing filter size for shift images
% OUTPUT
%   xyshift: Found shift in each Mask region in native (un-zoomed) pixels
%   xout, yout: Shifted Coordinates

rndsz = 64;
rndsz1 = 8;
xsz = ceil(max(x)/rndsz)*rndsz;
ysz = ceil(max(y)/rndsz)*rndsz;
zsz = ceil(max(z)/rndsz1)*rndsz1;

sz=[xsz ysz zsz];

cczoom = 4;
xsize = (sz(1)*cczoom);
ysize = (sz(2)*cczoom);
zsize = (sz(3)*cczoom);

x = single(cczoom*x);
y = single(cczoom*y);
z = single(cczoom*z);

mask0 = cormask==min(cormask);
%mask0 = cormask > 0;
refx = x(mask0);
refy = y(mask0);
refz = z(mask0);
if nargin>5
    [dx1,dy1,dz1] = findhistshift_3d(x,y,z,refx,refy,refz,xsize,ysize,zsize,sigma,cormask,stepsize*4);   
else
    [dx1,dy1,dz1] = findhistshift_3d(x,y,z,refx,refy,refz,xsize,ysize,zsize,sigma,cormask);
end
[xout1,yout1,zout1] = shiftcoords(x,y,z,cormask,[dx1,dy1,dz1]);

% convert back to Pixel Coordinates
xshift = dx1/cczoom;
yshift = dy1/cczoom;
zshift = dz1/cczoom;

xout = xout1/cczoom;
yout = yout1/cczoom;
zout = zout1/cczoom;

end

