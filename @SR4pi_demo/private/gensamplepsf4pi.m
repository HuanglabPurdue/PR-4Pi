function [samplepsf,startx,starty,startz,sizex,dz,dx] = gensamplepsf4pi(obj,psfsize,boxsize,phic,bin,Nzs,stagepos,chamberH,abertype)
prcal = obj.PSF4pi;
pixelsize = obj.Pixelsize;
zpos = linspace(-1.5,1.5,Nzs)';
xpos = zeros(Nzs,1);
ypos = zeros(Nzs,1);
I = ones(Nzs,4);
bg = zeros(Nzs,4);
x1 = cat(2,xpos,ypos,zpos,I,bg);
w = [1,1,1,1,1];
psfobj = prcal.psfobj;
PRstruct1 = prcal.PRstruct1;
PRstruct2 = prcal.PRstruct2;
psfobj.Pixelsize = pixelsize/bin;
psfobj.PSFsize = psfsize;
psfobj.Boxsize = boxsize*bin;
psfobj.Phi0 = phic;
if nargin > 7
    [psf4d_fit] = genpsf4pi_real(x1,w,psfobj,PRstruct1,PRstruct2,stagepos,chamberH,abertype);
else
    [psf4d_fit] = genpsf4pi_real(x1,w,psfobj,PRstruct1,PRstruct2);
end
psf_fit = [];
for ii = 1:4
    psf_fit = cat(2,psf_fit,squeeze(psf4d_fit(:,:,:,ii)));
end
samplepsf = cell(4,1);
for ii = 1:4
    f0 = psf4d_fit(:,:,:,ii);
    N = size(f0,1);
    Nz = size(f0,3);
    F = reshape(permute(f0,[2,1,3]),N*N*Nz,1);
    samplepsf{ii} = F;
end
startx = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
starty = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
startz = zpos(1);
sizex = size(psf4d_fit,1);
dz = zpos(2)-zpos(1);
dx = psfobj.Pixelsize;
end
