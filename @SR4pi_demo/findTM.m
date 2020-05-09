function findTM(obj)
% findTM - find affine transformations between the four quadrants.
% select a data for calculating the transformation matrix, data should
% be collected from the bottom emission path at low laser power 
% (we usually use an exposure time of 0.1 s). Sample can be either a
% flat cell structure or a few isolated beads taking up ~80% of the
% entire FOV. 
%
% Data format: 3D matrix (x,y,frame)
%
% Output:
%   Affs - affine transformation parameters, relative to the first quadrant
%          zm_all: zoom
%          trans_all: translation
%          ang_all: rotation
%          R: affine transformation matrix
%          invR: inverse affine transformation matrix
% 
% Output will be stored in obj and saved in file with key word '_FMTtransform' under obj.Resultpath\PSF measure\

%%    
resultpath = [obj.Resultpath,'\PSF measure\'];
offset = obj.Constantoffset;

[fileName, FileDir] = uigetfile([obj.Datapath,'*.mat'],'Select file:','MultiSelect','off');

%%
if ~exist(resultpath,'dir')
    mkdir(resultpath)
end
load(fullfile(FileDir,fileName))
h = dipshow(mean(cat(2,qd1,qd2,qd3,qd4),3));
diptruesize(h,200);
colormap(grey)
print(gcf,'-dpng','-r300',[resultpath, fileName(1:end-6),'sumprojection_image'])
close(h)
%% find Fourier-Mellin transform
zm_all = zeros(4,2);
trans_all = zeros(4,2);
ang_all = zeros(4,1);
R = zeros(3,3,4);
invR = zeros(3,3,4);
im1 = mean(qd1,3);
im1 = im1-offset;
im1(im1<0) = 1e-6;
im1 = im1./max(im1(:));
imref = repmat(im1,[1,4,1,1]);
imouts = [];
%dipshow(qdall);
for ii=1:4
    im2 = eval(['mean(qd', num2str(ii),',3);']);
    im2 = im2-offset;
    im2(im2<0) = 1e-6;
    im2 = im2./max(im2(:));
    rate = mean(im1(:))/mean(im2(:));
    im2 = im2*rate;

    [zm,trans,ang] = fmmatch(im2,im1);
    [out,R(:,:,ii)] = find_affine_trans(im2, im1, [[zm zm],trans,ang]);
    zm_fin=out(1:2);
    trans_fin=out(3:4);
    ang_fin=out(5);

    [imout] = double(affine_trans(im2,zm_fin,trans_fin,ang_fin));
    zm_all(ii,:)=zm_fin;
    trans_all(ii,:)=trans_fin;
    ang_all(ii)=ang_fin;
    [zm,trans,ang] = fmmatch(im1,im2);
    [~,invR(:,:,ii)] = find_affine_trans(im1, im2, [[zm zm],trans,ang]);
    imouts = cat(2,imouts,imout);
end
ov = joinchannels('RGB',imref,imouts);
dipshow(ov);
affS.zm_all = zm_all;
affS.trans_all = trans_all;
affS.ang_all = ang_all;
affS.R = R;
affS.invR = invR;
obj.Affs = affS;

save([resultpath, fileName(1:end-6), '_FMTtransform.mat'],'affS');
end



