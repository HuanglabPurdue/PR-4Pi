% for v14 iPALMast_analysisv14scmos.m

function [rms, rmp]=iPALMast_findmom_givenC_v1(subims,xf,yf,bgf,trans_dx,trans_dy)
N = size(subims,3);
Nqd = size(subims,4);
redmom = zeros(N,Nqd);
for ii=1:Nqd
    sigma=1.1;
    [x, y]=meshgrid(0:size(subims,1)-1,0:size(subims,1)-1);
    xxrep=repmat(x,[1 1 size(subims,3)]);
    yyrep=repmat(y,[1 1 size(subims,3)]);
    xfrep=repmat(reshape(xf+trans_dx(:,ii),[1 1 size(subims,3)]),[size(subims,1) size(subims,2)]);
    yfrep=repmat(reshape(yf+trans_dy(:,ii),[1 1 size(subims,3)]),[size(subims,1) size(subims,2)]);
    R=sqrt((xxrep-xfrep).^2+(yyrep-yfrep).^2);
    mask = exp(-R.^2./2./sigma./sigma).*R.^0;
    momarr = mask.*(subims(:,:,:,ii)-1.*bgf(:,:,:,ii)); % zeroth moment
    redmom(:,ii)=squeeze(sum(sum(momarr)));                                 % sum could be changed into a fitting with Gaussian
end
redmom(redmom<=0) = 1e-6;
rmp=(redmom(:,1)-redmom(:,3))./(redmom(:,1)+redmom(:,3));
rms=(redmom(:,4)-redmom(:,2))./(redmom(:,4)+redmom(:,2));