function [psf4dt_wT,psf4dt_noT,gainR4dt,cor1,trans_dx,trans_dy,rt_coeff] = makerois4pi_sCMOS(datapath,filename,obj,plotflag)
% select interference bead data folder
bxsz = obj.Boxsize;
qcenter = obj.Quadcenter;
gcal = load(obj.Gainpath);
thresh = obj.Peakthresh;
flipsigns=[0 0 1 1];
%qcenter =[388, 573, 1400, 1586];
[qd_offset] = iPALMscmos_makeqds(gcal.ccdoffset,qcenter,flipsigns);
[qd_gain] = iPALMscmos_makeqds(gcal.gain,qcenter,flipsigns);
[qd_var] = iPALMscmos_makeqds(gcal.ccdvar,qcenter,flipsigns);

load(fullfile(datapath,filename))
imsz = size(qd1);
qd_offsetL = repmat(qd_offset,[1,1,imsz(3),1]);
qd_gainL = repmat(qd_gain,[1,1,imsz(3),1]);
qd_gainR = squeeze(qd_var./qd_gain./qd_gain);

aveqd1 = (qd1-qd_offsetL(:,:,:,1))./qd_gainL(:,:,:,1);
aveqd2 = (qd2-qd_offsetL(:,:,:,2))./qd_gainL(:,:,:,2);
aveqd3 = (qd3-qd_offsetL(:,:,:,3))./qd_gainL(:,:,:,3);
aveqd4 = (qd4-qd_offsetL(:,:,:,4))./qd_gainL(:,:,:,4);    

[q1, q2, q3, q4]=align(aveqd1,aveqd2,aveqd3,aveqd4,obj.Affs);
sz = 16/5;    
rangemin=[17, 17];
rangemax=[150, 150];


sumim = q1+q4+q2+q3;
im_unif=unif(sumim,[sz sz 0],'rectangular')-unif(sumim,[2*sz 2*sz 0],'rectangular');
im_max=(im_unif>=.999*maxf(im_unif,[2*sz 2*sz 0],'rectangular'))&(im_unif>thresh);
centers = findcoord(im_max);
if plotflag == 1
    impeak = cHistRecon3D(imsz(1),imsz(2),imsz(3),single(centers(:,2)),single(centers(:,1)),single(centers(:,3)),0);
    Imax = max(sumim(:));
    overlay(sumim./Imax.*500,impeak)
end
    
cor1 = centers;
mask = cor1(:,1)<rangemax(1)&cor1(:,1)>rangemin(1)...
    &cor1(:,2)<rangemax(2)&cor1(:,2)>rangemin(2);
cor1 = cor1(mask,:);
Num = numel(cor1(:,1));

% find cor2, cor3 and cor4 from affine transform
zm = 1./obj.Affs.zm_all;
trans = obj.Affs.trans_all;
ang = obj.Affs.ang_all;

cc = imsz(1)/2;
Ta = [1 0 -cc-0.5;0 1 -cc-0.5;0 0 1];
Tb = [1 0 cc+0.5;0 1 cc+0.5;0 0 1];
rt = [1,0,0,1];
for nn = 2:4
    Rx = Tb*[cos(ang(nn)), sin(ang(nn)),0;-sin(ang(nn)),cos(ang(nn)),0;0,0,1]*Ta;
    Zx = Tb*[zm(nn,1),0,0;0,zm(nn,2),0;0,0,1]*Ta;
    Tx = [1,0,trans(nn,1);0,1,trans(nn,2);0,0,1];
    M = Tx*Zx*Rx;
    rt = cat(2,rt,M(1,1),M(1,2),M(2,1),M(2,2));
    tmp = M*[cor1(:,1),cor1(:,2),ones(Num,1)]';
    cor_tmp = [tmp(1,:)',tmp(2,:)',cor1(:,3)];
    eval(['cor',num2str(nn),' = cor_tmp;']);
end

rt_coeff = rt';
trans_dx = cat(2,cor1(:,1)-round(cor1(:,1)),cor2(:,1)-round(cor2(:,1)),cor3(:,1)-round(cor3(:,1)),cor4(:,1)-round(cor4(:,1)));
trans_dy = cat(2,cor1(:,2)-round(cor1(:,2)),cor2(:,2)-round(cor2(:,2)),cor3(:,2)-round(cor3(:,2)),cor4(:,2)-round(cor4(:,2)));

[qd1sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(q1));
[qd2sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(q2));
[qd3sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(q3));
[qd4sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(q4));

psf4dt_wT = cat(4,qd1sub,qd2sub,qd3sub,qd4sub);
if plotflag == 1
    im1 = cat(2,q3,q4);
    im2 = cat(2,q1,q2);
    colorim=joinchannels('RGB',im1,im2);
    h = dipshow(colorim);
    diptruesize(h,400)
end
clear q1 q2 q3 q4

[qd1sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(aveqd1));
[qd2sub]=cMakeSubregions(round(cor2(:,2))-1,round(cor2(:,1))-1,cor2(:,3),bxsz,single(aveqd2));
[qd3sub]=cMakeSubregions(round(cor3(:,2))-1,round(cor3(:,1))-1,cor3(:,3),bxsz,single(aveqd3));
[qd4sub]=cMakeSubregions(round(cor4(:,2))-1,round(cor4(:,1))-1,cor4(:,3),bxsz,single(aveqd4));

psf4dt_noT = cat(4,qd1sub,qd2sub,qd3sub,qd4sub);
clear aveqd1 aveqd2 aveqd3 aveqd4

qcL1 = repmat(qd_gainR(:,:,1),[1,1,imsz(3)]);
qcL2 = repmat(qd_gainR(:,:,2),[1,1,imsz(3)]);
qcL3 = repmat(qd_gainR(:,:,3),[1,1,imsz(3)]);
qcL4 = repmat(qd_gainR(:,:,4),[1,1,imsz(3)]);

[qc1sub]=cMakeSubregions(cor1(:,2)-1,cor1(:,1)-1,cor1(:,3),bxsz,single(qcL1));
[qc2sub]=cMakeSubregions(round(cor2(:,2))-1,round(cor2(:,1))-1,cor2(:,3),bxsz,single(qcL2));
[qc3sub]=cMakeSubregions(round(cor3(:,2))-1,round(cor3(:,1))-1,cor3(:,3),bxsz,single(qcL3));
[qc4sub]=cMakeSubregions(round(cor4(:,2))-1,round(cor4(:,1))-1,cor4(:,3),bxsz,single(qcL4));

gainR4dt = cat(4,qc1sub,qc2sub,qc3sub,qc4sub);

if plotflag == 1
    h = dipshow(cat(2,qd1sub,qd2sub,qd3sub,qd4sub));
    diptruesize(h,1000)
end


% im1 = zeros(imsz);
% im2 = zeros(imsz);
% im3 = zeros(imsz);
% 
% s1 = floor(bxsz/2);
% s2 = ceil(bxsz/2);
% cor2d = floor(cor2);
% 
% for ii = 1:Num
%     
%     im1(cor1(ii,2)-s1+1:cor1(ii,2)+s2,cor1(ii,1)-s1+1:cor1(ii,1)+s2,cor1(ii,3)+1) = psf4dt_noT(:,:,ii,1)+im1(cor1(ii,2)-s1+1:cor1(ii,2)+s2,cor1(ii,1)-s1+1:cor1(ii,1)+s2,cor1(ii,3)+1);
%     im2(cor2d(ii,2)-s1+1:cor2d(ii,2)+s2,cor2d(ii,1)-s1+1:cor2d(ii,1)+s2,cor2d(ii,3)+1) = psf4dt_noT(:,:,ii,2)+im2(cor2d(ii,2)-s1+1:cor2d(ii,2)+s2,cor2d(ii,1)-s1+1:cor2d(ii,1)+s2,cor2d(ii,3)+1);
%     
%     roi = double(affine_trans(psf4dt_noT(:,:,ii,2),[1,1],[trans_dx(ii,2),trans_dy(ii,2)],0));
%     im3(cor1(ii,2)-s1+1:cor1(ii,2)+s2,cor1(ii,1)-s1+1:cor1(ii,1)+s2,cor1(ii,3)+1) = roi+im3(cor1(ii,2)-s1+1:cor1(ii,2)+s2,cor1(ii,1)-s1+1:cor1(ii,1)+s2,cor1(ii,3)+1);
% end



