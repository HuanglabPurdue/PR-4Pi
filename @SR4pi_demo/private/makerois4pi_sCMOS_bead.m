function [psf4dt_wT,psf4dt_noT,gainR4dt,cor1,trans_dx,trans_dy,rot_coeff] = makerois4pi_sCMOS_bead(datapath,filename,obj,plotflag)
% select interference bead data folder
qcenter = obj.Quadcenter;
gcal = load(obj.Gainpath);
bxsz = obj.Boxsize;
flipsigns=[0 0 1 1];

[qd_offset] = iPALMscmos_makeqds(gcal.ccdoffset,qcenter,flipsigns);
[qd_gain] = iPALMscmos_makeqds(gcal.gain,qcenter,flipsigns);
[qd_var] = iPALMscmos_makeqds(gcal.ccdvar,qcenter,flipsigns);

tmp = load(fullfile(datapath,filename));
namei = fields(tmp);
sz = size(tmp.(namei{1}));
if numel(sz) == 3 
    sz(4) = 1;
end

for ii = 1:length(namei)
    qdi = reshape(tmp.(namei{ii}),[sz(1),sz(2),sz(3)*sz(4)]);
    eval(['qd',num2str(ii),'=qdi;'])
end
imsz = size(qd1);
qd_offsetL = repmat(qd_offset,[1,1,imsz(3),1]);
qd_gainL = repmat(qd_gain,[1,1,imsz(3),1]);
qd_gainR = squeeze(qd_var./qd_gain./qd_gain);

aveqd1 = (qd1-qd_offsetL(:,:,:,1))./qd_gainL(:,:,:,1);
aveqd2 = (qd2-qd_offsetL(:,:,:,2))./qd_gainL(:,:,:,2);
aveqd3 = (qd3-qd_offsetL(:,:,:,3))./qd_gainL(:,:,:,3);
aveqd4 = (qd4-qd_offsetL(:,:,:,4))./qd_gainL(:,:,:,4);    

[q1, q2, q3, q4]=align(aveqd1,aveqd2,aveqd3,aveqd4,obj.Affs);

qcL1 = repmat(qd_gainR(:,:,1),[1,1,imsz(3)]);
qcL2 = repmat(qd_gainR(:,:,2),[1,1,imsz(3)]);
qcL3 = repmat(qd_gainR(:,:,3),[1,1,imsz(3)]);
qcL4 = repmat(qd_gainR(:,:,4),[1,1,imsz(3)]);

sumim = q1+q4+q2+q3;
[cutim,centers] = pickbead(sumim,bxsz,1,'stack',[]);
    
cor1 = [centers(1).*ones(imsz(3),1),centers(2).*ones(imsz(3),1),[0:imsz(3)-1]'];
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

rot_coeff = rt';
trans_dx = cat(2,cor1(:,1)-round(cor1(:,1)),cor2(:,1)-round(cor2(:,1)),cor3(:,1)-round(cor3(:,1)),cor4(:,1)-round(cor4(:,1)));
trans_dy = cat(2,cor1(:,2)-round(cor1(:,2)),cor2(:,2)-round(cor2(:,2)),cor3(:,2)-round(cor3(:,2)),cor4(:,2)-round(cor4(:,2)));

qd1sub = pickbead(q1,bxsz,Num,'single',cor1);
qd2sub = pickbead(q2,bxsz,Num,'single',cor1);
qd3sub = pickbead(q3,bxsz,Num,'single',cor1);
qd4sub = pickbead(q4,bxsz,Num,'single',cor1);
psf4dt_wT = cat(4,qd1sub,qd2sub,qd3sub,qd4sub);

qd1sub = pickbead(aveqd1,bxsz,Num,'single',cor1);
qd2sub = pickbead(aveqd2,bxsz,Num,'single',round(cor2));
qd3sub = pickbead(aveqd3,bxsz,Num,'single',round(cor3));
qd4sub = pickbead(aveqd4,bxsz,Num,'single',round(cor4));
psf4dt_noT = cat(4,qd1sub,qd2sub,qd3sub,qd4sub);

qc1sub = pickbead(qcL1,bxsz,Num,'single',cor1);
qc2sub = pickbead(qcL2,bxsz,Num,'single',round(cor2));
qc3sub = pickbead(qcL3,bxsz,Num,'single',round(cor3));
qc4sub = pickbead(qcL4,bxsz,Num,'single',round(cor4));
gainR4dt = cat(4,qc1sub,qc2sub,qc3sub,qc4sub);


if plotflag == 1
    h = dipshow(cat(2,qd1sub,qd2sub,qd3sub,qd4sub));
    diptruesize(h,1000)
    im1 = cat(2,q3,q4);
    im2 = cat(2,q1,q2);
    colorim=joinchannels('RGB',im1,im2);
    h = dipshow(colorim);
    diptruesize(h,400)
end


