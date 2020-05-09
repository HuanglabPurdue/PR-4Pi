function [phic,p_phi,xp_sr,indm] = findphi0_loc(sr_data,phi_data,acc_thresh,obj,p_phi,plotflag)
srcenter = obj.Srcenter;
phi_offset = obj.Phi_offset;
R = 128;
sigma = 4;
[dataim,scale_data,xc_data,yc_data] = genblurimg(R,sr_data,phi_data,sigma);
dataim = dataim-mean(dataim(:));
dataim = dataim/std(dataim(:));
if isempty(p_phi)
    p_phi0 = obj.Astparamp.p_phi(1);
    %est = fminsearch(@(x)findslope4pi(x,p_phi0,dataim,R,sigma),-0.1,optimset('MaxIter',100,'Display','off'));
    [sse,refim,scale_ref,xc_ref,yc_ref,accx]=findslope4pi(0,p_phi0,dataim,R,sigma);
    L = double([1:numel(accx)]);
    [val,ind] = findpeaks(accx,L);
    vp = val(val>acc_thresh);
    xp = ind(val>acc_thresh);
    [xp1,vp1] = redefinepeak(xp,vp);
    xp_sr = (xp1+xc_data-(size(refim,1)-1)/2).*scale_data + min(sr_data);
%     tmp = cat(1,xp_sr,vp1);
%     tmp1 = sortrows(tmp',-2);
%    p_phi = 2*pi/mean(abs(tmp1(1,1)-tmp1(2,1)));
    maskvp = vp1>3000;
    if sum(maskvp)>1
        p_phi = -2*pi/mean(diff(xp_sr(maskvp)));
    else
        p_phi = p_phi0;
    end
else
    [sse,refim,scale_ref,xc_ref,yc_ref,accx]=findslope4pi(0,p_phi,dataim,R,sigma);
end

L = double([1:numel(accx)]);
[val,ind] = findpeaks(accx,L);
vp = val(val>acc_thresh);
xp = ind(val>acc_thresh);
if numel(xp)>1
    [xp1,vp1] = redefinepeak(xp,vp);
else
    xp1 = xp;
    vp1 = vp;
end
 
if plotflag == 1
    figure;plot(L,accx)
    hold on
    plot(xp1,vp1,'o')
end

xp_sr = (xp1+xc_data-(size(refim,1)-1)/2).*scale_data + min(sr_data);

vp_tmp = vp1;
[~,ind1] = max(vp_tmp);
% vp_tmp(ind1) = 0;
% [~,ind2] = max(vp_tmp);

phic = polyval([p_phi,0],srcenter)-p_phi*xp_sr(ind1);
% phic2 = polyval([p_phi,0],prcal.srcenter)-p_phi*xp_sr(ind2);
% w1 = abs(min([phic1+pi,pi-phic1,0]));
% w2 = abs(min([phic2+pi,pi-phic2,0]));
% 
% phic = (phic1*abs(w2) + phic2*abs(w1))/abs(w1-w2);
% indm = round((ind1*abs(w2) + ind2*abs(w1))/abs(w1-w2));
indm = ind1;
f1 = polyval([p_phi,0],sr_data)-p_phi*xp_sr(indm);
if plotflag == 1
    figure;plot(sr_data,phi_data,'.')
    hold on
    plot(sr_data,f1)
    plot(sr_data,f1+2*pi/3,'r')
    plot(sr_data,f1-2*pi/3,'r')
    plot(srcenter,phic,'o')
    ylim([-3,3])
end
phic = wrapToPi(phic+phi_offset);

function [xp1,vp1] = redefinepeak(xp,vp)
tmp = xp(1);
tmp1 = vp(1);
xp1 = [];
vp1 = [];
for ii = 1:numel(xp)-1
    if xp(ii+1)-xp(ii)<30
        tmp = cat(2,tmp,xp(ii+1));
        tmp1 = cat(2,tmp1,vp(ii+1));
        if ii+1 == numel(xp)
            xp1 = cat(2,xp1,mean(tmp));
            vp1 = cat(2,vp1,mean(tmp1));
        end
    else
        xp1 = cat(2,xp1,mean(tmp));
        vp1 = cat(2,vp1,mean(tmp1));
        tmp = xp(ii+1);
        tmp1 = vp(ii+1);
        if ii+1 == numel(xp)
            xp1 = cat(2,xp1,tmp);
            vp1 = cat(2,vp1,tmp1);
        end
    end
end


