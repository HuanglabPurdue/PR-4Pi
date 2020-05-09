% for v14 iPALMast_analysisv14scmos.m

function coords=iPALM_maxf_cand(imunif,sz,thresh)

im_max=(imunif>=.999*maxf(imunif,[sz sz 0],'rectangular'))&(imunif>thresh);
%im_max=(im_bgsub>=.999*maxf(im_bgsub,[sz sz 0],'rectangular'))&(im_bgsub>thresh);
%overim=overlay(imunif,im_max);
%dipshow(overim)
imsz=size(imunif,1);
a=find(im_max);
z=floor(a/imsz/imsz);
pnum=mod(a,imsz*imsz);
y=mod(pnum,imsz);
x=floor(pnum/imsz);
%        [out cds]=GPUFindCandidates(single(data(:,:,(st:ed))),PSFsigma,minI,ppt);  % PJC 5-4-10 correct bug if dimension 3 is 1
%         x=double(cds(1:end,1));
%         y=double(cds(1:end,2));
%         z=double(cds(1:end,3))+(jj-1)*NtimeBlock;
coords=[x y z];