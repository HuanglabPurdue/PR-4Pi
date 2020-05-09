function [sr_data,phi_data,cor_data,psf4dt_data,x0_data,z_phi,zast_data,sx_data,sy_data,maskall] = genini4pi_loc(psf4dt_wT,psf4dt_noT,cor1,ct_thresh,obj,trans_dx,trans_dy,p_phi,astfitflag,plotflag)
% I and bg
sz = size(psf4dt_noT);
Num = sz(3);
mask0 = ones(size(psf4dt_noT));
mask0(2:end-1,2:end-1,:,:) = 0;
tmp = psf4dt_noT;
tmp(mask0==0) = 0;
maski = mask0(:,:,1,:)==1;
L = sum(sum(mask0(:,:,1,1),1),2);
bg = zeros(1,1,Num,4);
for ii = 1:Num
    bgi = tmp(:,:,ii,:);
    bg1 = bgi(maski);
    bg2 = reshape(bg1,L,4);
    bg(:,:,ii,:) = median(bg2,1);
end
bg0 = squeeze(bg);
bgL = repmat(bg,[sz(1),sz(2),1,1]);
I = squeeze(sum(sum(sum(psf4dt_noT-bgL),1),2)).*1.5;
I0 = I;
I0(:,[1,3]) = repmat([I(:,1)+I(:,3)],1,2);
I0(:,[2,4]) = repmat([I(:,2)+I(:,4)],1,2);
I0(I0<100) = 100;
% x,y
datasum = squeeze(sum(psf4dt_wT,4));
datasum(datasum<=0) = 1e-6;
PSFsigma = 1.15;
iterations = 100;
fittype = 4;
[P, CRLB, LL] = GPUgaussMLEv2(single(datasum),PSFsigma,iterations,fittype);
x = P(:,2);
y = P(:,1);
[xx,yy] = meshgrid(0:sz(1)-1,0:sz(1)-1);
xxL = repmat(xx,[1,1,Num]);
yyL = repmat(yy,[1,1,Num]);
comx = squeeze(sum(sum(datasum.*xxL,1),2))./squeeze(sum(sum(datasum,1),2));
comy = squeeze(sum(sum(datasum.*yyL,1),2))./squeeze(sum(sum(datasum,1),2));
maskx = x>sz(1)-4|x<3;
masky = y>sz(1)-4|y<3;
x(maskx)=comx(maskx);
y(masky)=comy(masky);
% z localization
[rms, rmp] = iPALMast_findmom_givenC_v1(psf4dt_noT,x,y,bgL,trans_dx,trans_dy);
vec_amp = sqrt(rms.^2+rmp.^2);
rm1 = rms./vec_amp;
rm2 = rmp./vec_amp;
rm_complex = rm1-1i*rm2;
phi = angle(rm_complex);
z_phi = unwrap(phi).*obj.PhaseparamM0.zT./2./pi;
if plotflag == 1
    figure;
    scatter(rms,rmp,5,vec_amp)
    hold on
    %scatter(0,0,5,0)
    axis equal
end
sx = P(:,5);
sy = P(:,6);
sr = sx.^2 - sy.^2;
llr = -2*LL;

%% initial rejection
srlim = [-15,15];
maskc = vec_amp>ct_thresh & abs(rms)>1e-3 & abs(rms)<1-1e-3 & abs(rmp)>1e-3 & abs(rmp)<1-1e-3;
maskLL = llr<3e3;
if plotflag == 1
    figure;
    scatter(rms(maskc),rmp(maskc),5,vec_amp(maskc))
    hold on
    %scatter(0,0,5,0)
    axis equal
end
masksr = sr>srlim(1) & sr<srlim(2);

zfit0 = zeros(size(x));
start = sz(1)/2;
x0 = cat(2,x-start,y-start,zfit0,I0,bg0);

maskxy = abs(x0(:,1))<2.5 & abs(x0(:,2))<2.5;

if ~isempty(p_phi)
    N1 = numel(sr);
    density = zeros(N1,1);
    rphi = 0.2;
    parfor mm = 1:N1
        mask = ((sr-sr(mm)).^2+(phi-phi(mm)).^2./p_phi./p_phi) < rphi^2;
        density(mm) = sum(mask);
    end
    mask1 = density./N1./pi./rphi./rphi>0.005;
    maskall = maskxy & maskc & masksr & maskLL & mask1;
else
    maskall = maskxy & maskc & masksr & maskLL;
end


if ct_thresh == 0
    maskall = maskc>-1;
end

sr_data = sr(maskall);
phi_data = phi(maskall);
cor_data = cor1(maskall,:);
psf4dt_data = psf4dt_noT(:,:,maskall,:);
x0_data = x0(maskall,:);
sx_data = double(sx(maskall));
sy_data = double(sy(maskall));

if astfitflag == 1
    Nfit = numel(sr_data);
    zast_data = zeros(Nfit,1);
    astparam = obj.Astparamp;
    parfor ii=1:Nfit
        [zast_data(ii)] = iniz_astfit(sx_data(ii),sy_data(ii),astparam,0.1);
    end
end

if plotflag == 1
    figure;
    scatter(sx_data,sy_data,30,sr_data,'.')
    axis equal
    if astfitflag == 1
        figure;
        scatter(sx_data,sy_data,30,zast_data,'.')
        axis equal
    end
end



