function [mephi_ini,cpeak]=find_mephi_ini(dmap,srange,msz)

% normalization helps to emphasize on center map
normdmap=double(dmap.*gaussianblob(newim(size(dmap)),[size(dmap,1)/2 size(dmap,2)/2],size(dmap,1)/4,1));
[~,row,col]=findmax(normdmap);
cpeak=dmap(row,col);
[mephi_ini]=coords2mephi([col row],msz);
 
xx=linspace(-pi,pi,msz);
yy=linspace(-pi,pi,msz);
[xgrid, ygrid]=meshgrid(xx,yy);
rmap=sqrt((xgrid-mephi_ini(end,1)).^2+wrapToPi(ygrid-mephi_ini(end,2)).^2);
anglemap=angle((xgrid-mephi_ini(end,1))+(ygrid-mephi_ini(end,2))*1i);

maskr1=rmap<(srange(2))&rmap>(srange(1));

dv = 8;
angles = linspace(-pi,-pi/2,dv);
I = zeros(1,dv-1);
for ii = 1:dv-1
    anglemask = anglemap<angles(ii+1)&anglemap>angles(ii);
    mask = maskr1.*anglemask;
    patch = mask.*dmap;
    I(ii) = sum(patch(:));
end
[~,ind] = max(I);
da = 0.3;
angrange = [mean(angles(ind:ind+1))-da,mean(angles(ind:ind+1))+da];
anglemask = anglemap<angrange(2)&anglemap>angrange(1);
directmask=(xgrid-mephi_ini(end,1))<=0.01;
maskr=rmap<(srange(2)./3)&rmap>(srange(1)./3);
mask = maskr.*directmask.*anglemask;
patch = mask.*dmap;

%dipshow(patch);
[~,row,col]=findmax(patch);
mephi_ini(end+1,2)=yy(row);
mephi_ini(end,1)=xx(col);