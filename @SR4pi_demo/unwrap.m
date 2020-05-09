function unwrap(obj)
% unwrap - find z by unwrapping the found interference phase
%   It is based on a ridge finding algorithm developed by Huang in the
%   reference: Huang, F. et al. Ultra-High Resolution 3D Imaging of Whole Cells. Cell 166, 1028–1040 (2016)
%
% Output:
%   ResultCT.unwrap - localization result structre after unwrapping, besides the fields detailed in intifit, two additional fields are added
%       z - z positions, unit: nm
%       z_err - deviation from the ideal z position or the ridge line, unit: nm
%
% Output will be stored in obj and saved in a file with the keyword '_loc_z' under obj.Resultpath\CT fit\
%%
figfolder=[obj.Resultpath,'\CT fit\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end
%%
llrthresh = 800;
Ithresh = 700;
angctr = [-3, 3];% 3 sigma interval
bglim = [-3,4];
xLim = [1,obj.Imagesize].*obj.Pixelsize.*1e3;
nmed = obj.RI_med;
nimm = obj.RI_imm;
zthresh = [-800, 800]; % fit range from ast calibration
zasterr = [0, 1];
zT = obj.PhaseparamM0.zT*nimm/nmed;% modulation period, nm
zerrthresh = 80; %nm 
%%
res = obj.ResultCT.init;

mtc0 = (res.sx.^3./res.sy-res.sy.^3./res.sx);
normf = max(mean(mtc0)+1*std(mtc0),abs(mean(mtc0)-1*std(mtc0)));
mtc = mtc0./normf*pi/1;
res.mtc = mtc;
maskmtc = mtc<pi&mtc>-pi;% mtc mask
maskll = res.LL<llrthresh;% LL threshold
maskx = res.x>xLim(1)&res.x<xLim(2);

dv = 1e-3;
binc = [0.4:dv:2];
[angctnorm] = normalizedist(res.ct,binc,1);
[bgnorm] = normalizedist(res.bg,100,1);

maskzctr = angctnorm>angctr(1)&angctnorm<angctr(2);% contrast threshold
maskbg = bgnorm>bglim(1)&bgnorm<bglim(2);% background threshold
maskzast_err = res.zast_err<zasterr(2)&res.zast_err>zasterr(1);% astigmatism fitting threshold
maskzf = res.zast<zthresh(2)&res.zast>zthresh(1);% z range threshold
maskI = res.I>Ithresh;% photons range threshold
% apply
mask = maskI&maskzf&maskll&maskzast_err&maskzctr&maskbg&maskmtc&maskx;

freq = 2*pi/zT;
framenum = max(res.t)+1;

res1 = applymask(res,mask);


%% phase estimate for every optical section separately
phase_cyc = 1;

% get phase per cycle
fnum = phase_cyc.*framenum;
stopval = obj.Ridgepeak;
centermtc = [];
[currz, z_err, mephimask] = Mephi_z_4PiSMS(res1,fnum,stopval,centermtc,freq);

%% apply mask to every data
maskzerr = z_err<zerrthresh;% nm, adjust z resolution
maskall = maskzerr&(mephimask);
res1.z = -1.*currz;
res1.z_err = z_err;
res = applymask(res1,maskall);

obj.ResultCT.unwrap = res;
save([figfolder, obj.Savename, '_loc_z'],'res');

end