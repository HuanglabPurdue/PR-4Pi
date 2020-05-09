function rmoutlier(obj)
% rmoutlier - remove outliers from MLE localization results
%   a local density of each localization in a sr-z scatter plot was
%   calculated and localizations with a density less than obj.Mindensity
%   were removed. Ideally all localizations should fall into a single
%   stripe pattern in the sr-z scatter plot, multiple stripe patterns
%   indicates ghost localizations and those side stripes usually have a low
%   local density and can be therefore removed with a certain density
%   threshold.
%
% Output:
%   ResultPR.MLE.density - local density of each localization in a sr-z scatter plot
%   ResultPR.MLE.maskrm - filtering mask where outliers are set to zero
%
% Output will be stored in obj and saved in a file with the keyword '_mask' under obj.Resultpath\PR fit\
%%
figfolder=[obj.Resultpath,'\PR fit\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end
%%
res = obj.ResultPR.MLE.res;
sr = res.sx.^2-res.sy.^2;
z_fit = res.z;
Nr = numel(sr);

zT = median(obj.ResultPR.MLE.Phifit(:,1));
P_phi = obj.ResultPR.Phi0cali.P_phi;
rsr = 0.2;
if (Nr/max(res.stepN))<6e5
    dtype = 'smooth';
else
    dtype = 'discrete';
end
%%
if ~isfield(obj.ResultPR.MLE,'density')
    density = zeros(1,Nr);
    
    for nn = 1:max(res.stepN)
        
        maskstep = res.stepN == nn;
        z_fiti = z_fit(maskstep);
        sri = sr(maskstep);
        
        pz = median(P_phi)*zT/2/pi;
        N1 = numel(z_fiti);
        
        switch dtype
            case 'smooth'
                densityi = zeros(N1,1);
                
                tic
                parfor ii = 1:N1
                    mask = ((sri-sri(ii)).^2+(z_fiti-z_fiti(ii)).^2./pz./pz) < rsr^2;
                    densityi(ii) = sum(mask);
                end
                toc
            case 'discrete'
                densityi = zeros(N1,1);
                dz = 5;
                zrange = [-700:dz:700];
                srange = [-10:-dz/pz:10];
                [zz,ss] = meshgrid(zrange,srange);
                N2 = numel(zz);
                dcell = cell(N2,1);
                idcell = cell(N2,1);
                tic
                parfor ii = 1:N2
                    mask = sri>ss(ii)+dz/2/pz & sri<ss(ii)-dz/2/pz & z_fiti>zz(ii)-dz/2 & z_fiti<zz(ii)+dz/2;
                    % density(mask) = sum(mask);
                    dcell{ii} = ones(sum(mask),1).*sum(mask);
                    idcell{ii} = find(mask==1);
                end
                toc
                densityi(cell2mat(idcell)) = cell2mat(dcell);
        end
        
        density(maskstep) = densityi;
    end
    obj.ResultPR.MLE.density = density;
end
%%
density = obj.ResultPR.MLE.density;
maskall = zeros(1,Nr);
for nn = 1:max(res.stepN)
    maskstep = res.stepN == nn;
    N1 = sum(maskstep);
%     figure;ha = axes; scatter(sr(maskstep),z_fit(maskstep),1,density(maskstep),'.');
%     drawnow;
%     rect1 = getrect(ha);
%     zLim = [rect1(2),rect1(2)+rect1(4)];
    switch dtype
        case 'smooth'
        maskzsr = density./N1./pi./rsr./rsr>obj.Mindensity;
        maskzsr = maskzsr';
        case 'discrete'
        maskzsr = density./N1./dz./dz.*pz.*pz>obj.Mindensity;
        maskzsr = maskzsr;
    end
%     maskz = z_fit<zLim(2) & z_fit>zLim(1);
    maskzsr = maskzsr & maskstep;
    %figure;scatter(sr(maskzsr),z_fit(maskzsr),1,density(maskzsr),'.');
    maskall(maskzsr) = 1;
end
%%
maskzsr = maskall==1;
figure;scatter(sr(maskzsr),z_fit(maskzsr),1,density(maskzsr),'.');
obj.ResultPR.MLE.maskrm = maskzsr;

save([figfolder,obj.Savename,'_mask.mat'],'density','maskzsr','res')

end