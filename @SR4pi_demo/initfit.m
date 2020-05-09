function initfit(obj,dataname)
% initfit - initial localization for contrast algorithm, find x, y and interference phase
%   Data were collected at high laser power at 0.02 s exposure time, one
%   data file was saved for every 2000 frames, one dataset contains 10-90
%   data files. 
%
% Data format: each data file contains four 3D matrix with dimensions of
%              (x,y,frame) corresponding to data from the four imaging quandrants
%
% Output:
%   ResultCT.init.res - initial localization result structure
%           x - x positions, unit: nm
%           y - y positions, unit: nm
%           phi - interference phase
%           zast - z positions from astigmatism fitting with a Gaussian PSF model, unit: nm
%           I - photon of both objectives
%           bg - background per pixel from the sum of the four quadrants
%           LL - loglikelihood ratio
%           rms - modulation contrast from s-polarization
%           rmp - modulation contrast from p-polarization
%           ct -  modulation contrast from the central moment of the 4PiPSF patterns, ct = sqrt(rms^2+rmp^2)
%           sx - width of the astigmatism PSFs along the x-axis
%           sy - width of the astigmatism PSFs along the y-axis
%           zast_err - goodness of astigmatism fitting of the axial positions using a Gaussian PSF model, localizations with zast_err>1 will be removed during unwrap()
%           f - data file number
%           t - frame number relative to each data file
%           stepN - optical section number
%           cycN - cycle number of imaging through one cycle of all optical sections, it is the same as f for single-optical section imaging
%
% Output will be stored in obj and saved in a file with the keyword '_loc_ini_phi' under obj.Resultpath\CT fit\
%%
figfolder=[obj.Resultpath,'\CT fit\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end
%%
datafile = dir([obj.Datapath,dataname,'*.mat']);

%%
[res] = loc4pi_ct(obj.Datapath,datafile,figfolder,obj);
obj.Savename = datafile(1).name(1:end-4);
res.zast = res.zast - 1000;
obj.ResultCT.init = res;
save([figfolder, datafile(1).name(1:end-4), '_loc_ini_phi'],'res');

end