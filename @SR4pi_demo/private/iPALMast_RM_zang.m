% for v14 iPALMast_analysisv14scmos.m

% Version numbers:
% July30th 2014: now crop each image and their noisemap accordingly in
% their native scmos qd images.

% previous versions could be located from google drive versions.

function [z_ang, rms, rmp,ang_contrast, subims, z_angd]=iPALMast_RM_zang(obj,q1,q2,q3,q4,subsz,varim,gainim,res)
xf = res.x;
yf = res.y;
tlz_f = res.cor;
bg_f = res.bg;
affS = obj.Affs;
phaseparam = obj.PhaseparamM0;
affS.R(:,:,1) = eye(3);
varim = repmat(permute(varim,[1 2 4 3]),[1 1 size(q1,3) 1]);
gainim = repmat(permute(gainim,[1 2 4 3]),[1 1 size(q1,3) 1]);

xc = double(xf+tlz_f(:,2));
yc = double(yf+tlz_f(:,1));
tc = double(tlz_f(:,3));
N = numel(tc);
subims = zeros(subsz,subsz,N,4);
xcenter = zeros(N,4);
ycenter = zeros(N,4);
subvar = zeros(subsz,subsz,N,4);
subgain = zeros(subsz,subsz,N,4);

for ii=1:4
    eval(['im=q' num2str(ii) ';']);
    cocurr = [xc' ; yc'; ones(1,numel(xc))];
    transco = affS.R(:,:,ii)*cocurr;
    trans_x = transco(1,:);
    trans_y = transco(2,:);
    [subims(:,:,:,ii), t, l] = cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(im,[1 2 3])));
    xcenter(:,ii) = trans_x(:)-l;
    ycenter(:,ii) = trans_y(:)-t;
    
    [subvar(:,:,:,ii)] = cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(varim(:,:,:,ii),[1 2 3])));  % cmake subregions should crop noise maps
    [subgain(:,:,:,ii)] = cMakeSubregions(round(trans_y(:)),round(trans_x(:)),tc,subsz,single(permute(gainim(:,:,:,ii),[1 2 3])));  % cmake subregions should crop noise maps
end

sigma = 0.9;
sz = size(subims,1);
N = numel(xf);
models = zeros(sz,sz,N,4);
maxN = 5000;
repN = floor(N/maxN);
vec = cat(1,reshape(repmat([1:repN],maxN,1),repN*maxN,1),(repN+1).*ones(N-repN*maxN,1));
for ii = 1:4
    model = [];
    for nn = 1:repN+1
        maski = vec==nn;
        sampN = sum(maski);
        modeli=GPUgenerateBlobs(single(sz*sampN),single(sz),single(ycenter(maski,ii)+[0:sampN-1]'.*sz),single(xcenter(maski,ii)),single(ones(sampN,1)),single(ones(sampN,1).*sigma),single(ones(sampN,1).*sigma),single(zeros(sampN,1)),1);
        tmp = reshape(modeli',sz,sz,sampN);
        tmp = permute(tmp,[2,1,3]);
        model = cat(3,model,tmp);
    end
models(:,:,:,ii) = model;
end

bgs = reshape(repmat(bg_f',sz*sz,1),sz,sz,N);
int_est = zeros(N,4);
var_all=subims+subvar./subgain./subgain;

for ii = 1:4
    nom=(subims(:,:,:,ii)-bgs./4).*models(:,:,:,ii)./var_all(:,:,:,ii);
    denom=models(:,:,:,ii).^2./var_all(:,:,:,ii);
    int_est(:,ii)=squeeze(sum(sum(nom,1),2))./squeeze(sum(sum(denom,1),2))./2;% devided by 2 is possibly not necessary
end

rms=(int_est(:,1)-int_est(:,3))./(int_est(:,1)+int_est(:,3));% normalize to the photon count of s polarization
rmp=(int_est(:,4)-int_est(:,2))./(int_est(:,4)+int_est(:,2));% normalize to the photon count of p polarization
rm1=rms./sqrt(rms.^2+rmp.^2);
rm2=rmp./sqrt(rms.^2+rmp.^2);

phi_s = phaseparam.phis;
phi_p = phaseparam.phip;
c=(rm1-1i*rm2);                                                            
z_angd=angle(c);
ang2=[];
parfor ii=1:numel(rms)
    [ang2(ii,:)]=iPALM_est_angle(double(rms(ii)),double(rmp(ii)),phi_s,phi_p,wrapToPi(double(z_angd(ii)-phi_s)));
end

z_ang=wrapToPi(ang2(:,2));
ang_contrast=ang2(:,1);




