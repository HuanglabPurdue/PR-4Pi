%------ example code for using SR4pi_demo class on analyzing 4Pi-SMSN data------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Note: this example code demonstrates the usage of SR4pi_demo class, 
%       type 'help SR4pi_demo' in the command window for more information
%
% References: 
%
% (C) Copyright 2020                Huang Lab, Weldon School of Biomedical Engineering, 
%     All rights reserved           Purdue University, West Lafayette, IN, USA
%
%                                   
% Author: Sheng Liu, March 2020
%%
clearvars
addpath('PSF Toolbox\');
addpath('mex\')
addpath('test data\')
%% setup data path and result path
datapath = 'test data\';                                % directory of the data to be analyzed
resultpath = 'output\';                                 % directory of localization result
if ~exist(resultpath,'dir')
    mkdir(resultpath)
end
%% prepare parameters for localization, run time is ~5 minutes
srobj = SR4pi_demo(datapath,resultpath);                % create object from SR4pi_demo class
srobj.Gainpath = 'test data\gainmap.mat';               % directory of the gain calibration file
srobj.NA_bot = 1.4;                                     % numerical aperture of the bottom objective used for phase retrieval
srobj.Lambda = 0.675;                                   % emission wave length used for 4PiPSF model, unit: micron
srobj.Lambda_bead = 0.675;                              % emission wave length used for bead data, it is the wavelength relative to air and will be a constant after optimizing the hyper parameters (see optimparam), unit: micron
srobj.RI_imm = 1.516;                                   % refractive index of the immersion medium
srobj.RI_med = 1.351;                                   % refractive index of the sample medium
srobj.Pixelsize = 0.129;                                % pixel size relative to the sample plane, unit: micron
srobj.Peakthresh = 20;                                  % intensity peak of candidate emitters
% user input required for findTM(): 
%   select file: test data\cell_FMT_000_000.mat
srobj.findTM()                                          % find affine transformations between the four quadrants.
% user input required for caliparam(): 
%   select file: test data\bead_4pi_000_000.mat
%   select the emitter center from the popup window
srobj.caliparam()                                       % calibrate parameters of incoherent (astigmatism) and interferometric PSFs.
% user input required for getPR('bottom'): 
%   select file: test data\bead_bot_000_020.mat 
%   select the emitter center from the popup window
srobj.getPR('bottom')                                   % generate phase retrieved pupil function of the bottom emission path
% user input required for getPR('top'): 
%   select file: test data\bead_top_000_020.mat 
%   select the emitter center from the popup window
srobj.NA_top = 1.28;                                    % numerical aperture of the top objective used for phase retrieval,this value is sample dependent, it is normally less than 1.35, the value is selected so that the edge of the pupil function is smooth
srobj.getPR('top')                                      % generate phase retrieved pupil function of the top emission path
srobj.buildPSF()                                        % generate PR-4PiPSF model from phase retrieved pupil functions
% user input required for optimparam():
%   only for the first time calling optimparam()
%   select file: test data\bead_4pi_000_000.mat, 
%   select the emitter center from the popup window
srobj.optimparam('evaluate')                            % evaluate the goodness of initial hyper parameters
srobj.optimparam('lambda')                              % optimize emission wavelength (relative to air) for bead data by minimizing axial localization deviation within a defined z range
srobj.optimparam('Iratio')                              % optimize intensity ratio between top and bottom emission path by matching the slope of sr-z curves from data and PSF models
srobj.optimparam('srcenter')                            % optimize the sr (shape metric) for infocus PSF by minimizing the offset of sr-z curves from data and PSF models
srobj.optimparam('evaluate')                            % evaluate the goodness of optimized hyper parameters
% save srobj before localization
save([srobj.Resultpath,'tom20_srobj.mat'],'srobj')

%% test PR-4pi fitting, run time is ~5 minutes
dataname = 'cell1';                                     % keyword of data files
srobj.caliphi0(dataname,'single')                       % calibrate cavity phase for dataset from single-optical section imaging
srobj.mlefit(dataname)                                  % localization using maximum likelihood estimation (MLE) based on PR-4PiPSF model
srobj.Mindensity = 0.04;                                % minimum density threshold used in outlier removal on a sr-z density map, as the test data size is small, we used a higher threshold here, normally it is ~0.01.
srobj.rmoutlier()                                       % remove outliers from MLE localization results 
srobj.driftcorrection('PR')                             % post drift correction of localization results from PR-4Pi fitting

%% test contrast fitting, run time is ~2 minutes
dataname = 'cell1';                                     % keyword of data files
srobj.initfit(dataname);                                % initial localization, find x, y and interference phase
srobj.Ridgepeak = 0.1;                                  % stop value for ridge finding algorithm in contrast method
srobj.unwrap();                                         % find z by unwrapping the found interference phase
srobj.driftcorrection('CT')                             % post drift correction of localization results from contrast fitting

%% generate high resolution image of the 2D projection view of the whole data, color coded in z
srobj.Renderparam.zoomf = 10;                           % zoom factor of the high resolution image
srobj.Renderparam.pmax = 1;                             % maximum pixel value used for adjusting the image contrast
srobj.Renderparam.sigma = 1.5;                          % width of the Gaussian blur of each localization point, it is an average localization precision, unit: pixelsize/zoomf
% render result from PR-4Pi fitting
srobj.filtering('PR')                                   % filter out localizations with low photon count, high background and high LLR
srobj.renderSR('PR')                                    % generate and save high resolution image
% render result from contrast fitting
srobj.filtering('CT')                                   % filter out localizations with low photon count, high background and high LLR
srobj.renderSR('CT')                                    % generate and save high resolution image

%% output final result structure after drift correction
resp = srobj.ResultPR.Final.res;                        % from PR-4Pi fitting
resc = srobj.ResultCT.Final.res;                        % from contrast fitting
%% output result after final filtering step, which is used for rendering. 
resp_f = srobj.ResultPR.Final.filtered;                 % from PR-4Pi fitting
resc_f = srobj.ResultCT.Final.filtered;                 % from contrast fitting

