function renderSR(obj,type)
% renderSR - generate high resolution image of the 2D projection view of the whole data
%
% Input:
%   type - 'PR': for localizations from PR-4Pi algorithm
%          'CT': for localizations from contrast algorithm
%
% Generated high resolution image will be saved under obj.Resultpath as a .png file
%%
zoomi = obj.Renderparam.zoomf;
pmax = obj.Renderparam.pmax;
sigma = obj.Renderparam.sigma;
pxsz = obj.Pixelsize.*1e3;
xsz = ceil(max([diff(obj.Limit.x),diff(obj.Limit.y)])/pxsz);

res_sub = obj.(['Result',type]).Final.filtered;

sz1 = floor(diff(obj.Limit.x)/pxsz*zoomi);
sz2 = floor(diff(obj.Limit.y)/pxsz*zoomi);
resc.x = res_sub.x - obj.Limit.x(1);
resc.y = res_sub.y - obj.Limit.y(1);
resc.z = res_sub.z;
resc.z([1,2]) = obj.Limit.z;
[render_img] = renderslice1(resc,xsz,zoomi,pxsz,pmax,sigma,'jet');
img1 = render_img(1:sz2,1:sz1,:);


tmp = img1(img1>2);
[val,ed] = histcounts(tmp,100,'Normalization','cdf');
indct = find(val>0.99);
normf = ed(indct(1));

h = figure('Position',[100,100,1000,1000]);
h.InvertHardcopy = 'off';
ha = axes('Position',[0,0,1,1]);
image(img1./normf);
axis equal 
axis off
ha.YDir = 'normal';

imwrite(img1./normf,[obj.Resultpath,obj.Savename,'_',type,'.png'])


end