function [xshift,yshift,zshift] = findhistshift_3d(x,y,z,refx,refy,refz,xsize,ysize,zsize,sigma,cormask,stepsize)

refim = cHistRecon3D(xsize,ysize,zsize,refx,refy,refz,0);
refim = gaussf(permute(dip_image(stretch(refim),'uint8'),[2 1 3]),sigma);
ind = [min(cormask):max(cormask)];
ns = numel(ind);
ims = newim(xsize,ysize,zsize,ns);
xshift = zeros(ns,1);
yshift = zeros(ns,1);
zshift = zeros(ns,1);
for nn = 1:ns
    xi = single(x(cormask==ind(nn)));
    yi = single(y(cormask==ind(nn)));
    zi = single(z(cormask==ind(nn)));
    tempim = cHistRecon3D(xsize,ysize,zsize,xi,yi,zi,0);
    tempim = gaussf(permute(dip_image(stretch(tempim),'uint8'),[2 1 3]),sigma);
    ims(:,:,:,nn-1) = tempim;
    if nargin>11
       [shift1, shift2, shift3]=drift_correction_core3d_stack(refim,tempim,[0,0,round(stepsize*(nn-1))]);
    else
        [shift1, shift2, shift3]=drift_correction_core3d(refim,tempim,[0,0,0]);
    end
%     dxyz = findshift(refim,tempim,'iter');
%     xshift(nn) = dxyz(1);
%     yshift(nn) = dxyz(2);
%     zshift(nn) = dxyz(3);
    xshift(nn) = shift1;
    yshift(nn) = shift2;
    zshift(nn) = shift3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
end

xshift = xshift - xshift(1);
yshift = yshift - yshift(1);
zshift = zshift - zshift(1);

end