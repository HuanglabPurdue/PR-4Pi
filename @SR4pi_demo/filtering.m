function filtering(obj,type)
% filtering - filter out localizations with low photon count, high background and high LLR
%       Photon threshold was calculated at the peak of the photon distribution 
%   from all localizations and localizations with a photon count below the 
%   threshold were removed. 
%       Background threshold was estimated at mu+3sigma from 
%   a Gaussian fitting (with mean at mu and standard deviation at sigma) of the 
%   background distribution from all localizations and localizations with a 
%   background above the threshold were removed. 
%       LLR threshold was estimated at the 9% of the peak value of the LLR
%   distribution and localizations with a LLR above the threshold were removed.
%       z range was set at the 1% and 99.5% of the cumulative distribution of the
%   z localizations. Localizations outside the z range were removed
%
% Input:
%   type - 'PR': for localizations from PR-4Pi algorithm
%          'CT': for localizations from contrast algorithm
%
% Output:
%   Limit - threshold used in the filtering step, including threshold in x,y,z,I,bg,LL
%
% Output will be stored in obj
%%
res = obj.(['Result',type]).Final.res;
switch type
    case 'PR'
        res.I = res.I.*2;
        res.bg = res.bg.*4;
        llrange = [1000:20:4000];
    case 'CT'
        llrange = [10:5:800];
end

%% set photon limit
h = figure;
ha = axes;hold on
Imax = 10000;

hs = histogram(res.I,linspace(0,Imax,70));hs.Visible = 'off';
binc = (hs.BinEdges(1:end-1)+hs.BinEdges(2:end))/2;
tmp = sortrows(cat(1,hs.Values,binc)',-1);

[Ilim,idm] = min(tmp(1:2,2));
countm = tmp(idm,1);
hp = plot(binc,hs.Values,'parent',ha);

plot(Ilim,countm,'o','color',hp.Color,'parent',ha)
axis(ha,'tight')
ha.XLabel.String = 'photon';
ha.YLabel.String = 'count';
obj.Limit.I = Ilim;

%% set bg limit
h = figure;
ha = axes;hold on
bgmax = 200;

[~,mu,sigma] = normalizedist(res.bg,linspace(0,bgmax,100),0);
hs = histogram(res.bg,linspace(0,bgmax,100));hs.Visible = 'off';
binc = (hs.BinEdges(1:end-1)+hs.BinEdges(2:end))/2;
bglim = mu+3*sigma;
[~,idm] = min(abs(bglim-binc));
countm = hs.Values(idm);
hp = plot(binc,hs.Values,'parent',ha);
plot(bglim,countm,'o','color',hp.Color,'parent',ha)
axis(ha,'tight')
ha.XLabel.String = 'background';
ha.YLabel.String = 'count';
obj.Limit.bg = bglim;

%% set LLR limit
h = figure;
ha = axes;hold on

llc = 0.09;

[N,ed] = histcounts(res.LL,llrange,'Normalization','probability');
id = find((N./max(N))>=llc);
llim = (ed(id(end)) + ed(id(end)+1))/2;

hs = histogram(res.LL,llrange);hs.Visible = 'off';
binc = (hs.BinEdges(1:end-1)+hs.BinEdges(2:end))/2;

[~,idm] = min(abs(llim-binc));
countm = hs.Values(idm);
hp = plot(binc,hs.Values,'parent',ha);
plot(llim,countm,'o','color',hp.Color,'parent',ha)
axis(ha,'tight')
ha.XLabel.String = 'LLR';
ha.YLabel.String = 'count';
obj.Limit.LL = llim;

%% set z range
zlimi = zeros(1,2);
pz = 0.005;
[N,ed] = histcounts(res.z,100,'Normalization','cdf');

id = find(N<=pz);
zlimi(1) = (ed(id(end)) + ed(id(end)+1))/2;

id = find(N>=1-pz);
zlimi(2) = (ed(id(1)) + ed(id(2)))/2;

h = figure;
ha = axes;hold on
hs = histogram(res.z,100);hs.Visible = 'off';
binc = (hs.BinEdges(1:end-1)+hs.BinEdges(2:end))/2;
idm = zeros(1,2);
[~,idm(1)] = min(abs(zlimi(1)-binc));
[~,idm(2)] = min(abs(zlimi(2)-binc));

countm = hs.Values(idm);
hp = plot(binc,hs.Values,'parent',ha);
plot(zlimi,countm,'o','color',hp.Color,'parent',ha)
axis(ha,'tight')
ha.XLabel.String = 'z (nm)';
ha.YLabel.String = 'count';

obj.Limit.z = zlimi;

%% set xy range for full FOV
obj.Limit.x = [1,obj.Imagesize].*obj.Pixelsize.*1e3;      % nm
obj.Limit.y = [1,obj.Imagesize].*obj.Pixelsize.*1e3;      % nm

%% apply filter 

[masksub] = genmask(res,obj.Limit);
res_sub = applymask(res,masksub);

obj.(['Result',type]).Final.filtered = res_sub;
end






