% for v19 iPALMast_analysisv14scmos.m
% v19 
% includes sCMOS noise

function [subims, coords, locmaxc, subvar_g]=iPALMast_sumim_seg(im,thresh,subsz,varim,sumvg,thresholdflag)

imunf=iPALM_unif_sCMOS(im,varim,5);
% imunf=iPALM_unif_EMCCD(im,5);
maxfsz=5;
locmaxc=iPALM_maxf_cand(unif(imunf,[2 2 0]),maxfsz,thresh);
%impeak = cHistRecon3D(imsz(1),imsz(2),imsz(3),single(locmaxc(:,2)),single(locmaxc(:,1)),single(locmaxc(:,3)),0);
%overlay(im,impeak)
rangemin=[17 17];
rangemax=[150 150];
if thresholdflag==1
    mask=locmaxc(:,1)<rangemax(1)&locmaxc(:,1)>rangemin(1)...
        &locmaxc(:,2)<rangemax(2)&locmaxc(:,2)>rangemin(2);
    locmaxc=locmaxc(mask,:);
end
% subsz=15;
if isempty(locmaxc)
    subims=[];
    coords=[];
else
    
[subims, x, y]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(permute(im,[1 2 3])));
[subvar_g]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(repmat(permute(sumvg,[1 2 3]),[1 1 size(im,3)])));
% [subg]=cMakeSubregions(locmaxc(:,2),locmaxc(:,1),locmaxc(:,3),subsz,single(repmat(permute(gainim,[1 2 3]),[1 1 size(im,3)])));

coords = [x y locmaxc(:,3)];
end