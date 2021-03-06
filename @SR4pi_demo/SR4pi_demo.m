classdef SR4pi_demo < handle
    % (C) Copyright 2020                Huang Lab, Weldon School of Biomedical Engineering, 
    %     All rights reserved           Purdue University, West Lafayette, IN, USA
    %
    %                                   
    % Author: Sheng Liu, March 2020
    %
    % SR4pi_demo class for localization of 4Pi data from 4Pi-SMSN
    % microscope systems. Two localization methods are provided:
    % Method1: PR-4Pi algorithm
    %          Based on PR-4PiPSF models generated from coherent pupil
    %          functions obtained from phase retrieval.
    %          Reference: Liu, S. et al. Enhanced 4Pi single-molecule localization microscopy with coherent pupil based localization. Communications Biology, (2020) 
    %
    % Method2: contrast algorithm
    %          Axial localization is based on the modulation contrast of
    %          the central moment of 4PiPSF patterns and lateral localization 
    %          is based on incoherent PSF patterns.
    %          Reference: Huang, F. et al. Ultra-High Resolution 3D Imaging of Whole Cells. Cell 166, 1028�1040 (2016)
    %
    % Create object: srobj = SR4pi_demo(datapath,resultpath);
    % 
    % SR4pi_demo Methods:
    %   findTM - find affine transformations between the four quadrants.
    %   caliparam - calibrate parameters of incoherent (astigmatism) and interferometric PSFs.
    %   getPR - generate phase retrieved pupil functions
    %   buildPSF - generate PR-4PiPSF model from phase retrieved pupil functions
    %   optimparam - optimize hyper parameters for PR-4PiPSF model, including lambda, Iratio and srcenter
    %   driftcorrection - post drift correction of localization results
    %   filtering - filter out localizations with low photon count, high background and high LLR
    %   renderSR - generate high resolution image of the 2D projection view of the whole data 
    %
    %   Methods for PR-4Pi algorithm:
    %   caliphi0 - calibrate cavity phase, which is assciated with room temperature and humidity change
    %   mlefit - localization using maximum likelihood estimation (MLE) based on PR-4PiPSF model
    %   rmoutlier - remove outliers from MLE localization results 
    %
    %   Methods for contrast algorithm:
    %   initfit - initial localization, find x, y and interference phase
    %   unwrap - find z by unwrapping the found interference phase
    %
    % Naming convention:
    %   Properties: start with a uppercase letter (field name of a struct property uses either uppercase or lowercase letter)
    %   Methods: start with a lowercase letter
    %
    properties
        Constantoffset = 100;                       % a constant offset of the raw camera frame, used only for calibration
        Constantgain = 2;                           % a constant gain of the raw camera frame, used only for calibration
        % affine transformation parameters, relative to the first quadrant
            % zm_all: zoom
            % trans_all: translation
            % ang_all: rotation
            % R: affine transformation matrix
            % invR: inverse affine transformation matrix
        Affs;
        Quadcenter = [282,490,1673,1881];           % horizontal center of each quadrant relative to the full camera frame
        Boxsize = 16;                               % subregion size of selected emitters for localization
        Resultpath;                                 % directory of localization result
        Savename;                                   % key word of saved localization result files
        Datapath;                                   % directory of the data to be analyzed
        Gainpath;                                   % directory of the gain calibration file
        NA_bot = 1.4;                               % numerical aperture of the bottom objective used for phase retrieval
        NA_top = 1.35;                              % numerical aperture of the top objective used for phase retrieval
        Lambda = 0.675;                             % emission wave length used for 4PiPSF model, unit: micron     
        Lambda_bead = 0.675;                        % emission wave length used for bead data, it is the wavelength relative to air and will be a constant after optimizing the hyper parameters (see optimparam), unit: micron
        RI_imm = 1.516;                             % refractive index of the immersion medium
        RI_med = 1.351;                             % refractive index of the sample medium
        Pixelsize = 0.129;                          % pixel size relative to the sample plane, unit: micron
        Srcenter;                                   % shape metric for infocus astigmatism PSF
        Phi_offset;                                 % offset of measured cavity phase to the theoretical cavity phase, it is near zero, unit: radian
        Iratio = 0.6;                               % transmission ratio between top and bottom emission path, from 0 to 1 
        Modulationdepth = 0.7;                      % modulation strength of interferometric PSFs,from 0 to 1 
        OTFsigma = 1.6;                             % sigma of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaX), unit: micron
        % calibration parameters for Gaussian fitting of astigmatism PSFs, used in contrast method
            % estx: coefficients of the calibration curve along the x-axis
            % esty: coefficients of the calibration curve along the y-axis
            % zstep: axis positions used for calibration
            % sx: width of the astigmatism PSFs along the x-axis
            % sy: width of the astigmatism PSFs along the y-axis
        Astparamc;    
        % calibration parameters for Gaussian fitting of astigmatism PSFs, used in PR-4Pi method
            % estx: coefficients of the calibration curve along the x-axis
            % esty: coefficients of the calibration curve along the y-axis
            % sx: width of the astigmatism PSFs along the x-axis
            % sy: width of the astigmatism PSFs along the y-axis
            % sr: shape metric, sr=sx^2-sy^2
            % p_phi: coefficients from linear fitting of the sr-phi curve, phi is the interference phase
            % p: coefficients from linear fitting of the sr-z curve, z is the axial positions
        Astparamp;    
        % calibration parameters of the central moment modulation of the 4PiPSFs
            % phis: relative phase in s-polarization, unit: radian
            % phip: relative phase in p-polarization, unit: radian
            % phi: interference phase, unit: radian
            % zT: modulation period, unit: nm
            % ang_coeff: coefficients of fitting the modulation curve with a cos function, zT = 2pi/ang_coeff(1)
        PhaseparamM0;
        % calibration parameters of the third moment modulation of the 4PiPSFs
            % phis: relative phase in s-polarization, unit: radian
            % phip: relative phase in p-polarization, unit: radian
            % phi: interference phase, unit: radian
            % zT: modulation period, unit: nm
            % ang_coeff: coefficients of fitting the modulation curve with a cos function, zT = 2pi/ang_coeff(1)
        PhaseparamM3;
        PSFcalidata;                                % 4PiPSF data used for calibration, including a z-stack PSFs for each quadrant, see caliparam
        PSFoptidata;                                % 4PiPSF data used for optimizing the hyper parameters
        PSFtop;                                     % object of OptimPR_Ast class for phase retrieval of top emission path
        PSFbot;                                     % object of OptimPR_Ast class for phase retrieval of bottom emission path
        % 4PiPSF model
            % psfobj: object of PSF_4pi class
            % PRstruct1: major parameters for top emission path
            % PRstruct2: major parameters for bottom emission path
        PSF4pi;    
        PSFlib;                                     % 4PiPSF library used for localization
        Maxframe = 2000;                            % frame number per file
        Peakthresh = 20;                            % intensity peak of candidate emitters
        % localization results from PR-4Pi algorithm
            % Phi0cali: calibration results of cavity phase
            % MLE: localization results from maximum likelihood estimation 
            % Final: localization results after outlier removal and drift correction
        ResultPR;
        % localization results from contrast algorithm 
            % init: localization results of x,y and interference phase
            % unwrap: localization results after unwrapping the interference phase
            % Final: localization results after drift correction
        ResultCT;    
        IterationNum = 40;                          % number of iterations used in MLE
        Dampingfactor = 1;                          % damping factor for Levenberg-Marquardt algorithm used in MLE
        Stepsize = 0.6;                             % scanning step size used in multi section imaging, unit: micron
        Imagesize = 168;                            % image size of each quadrant
        Ridgepeak = 0.15;                           % stop value for ridge finding algorithm in contrast method
        Paramcheck;                                 % goodness of optimized parameters: srcenter, lambda, Iratio
        Mindensity = 0.01;                          % minimum density threshold used in outlier removal on a sr-z density map
        % parameters used for generating high resolution image
            % zoomf: zoom factor of the high resolution image
            % pmax: maximum pixel value used for adjusting the image contrast
            % sigma: width of the Gaussian blur of each localization point, it is an average localization precision, unit: pixelsize/zoomf
        Renderparam;
        Limit;                                      % threshold used in the filtering step, including threshold in x,y,z,I,bg,LL
    end
    
    methods
        
        function obj = SR4pi_demo(datapath,resultpath)
            obj.Datapath = datapath;
            obj.Resultpath = resultpath;
        end
        
        function set.Lambda(obj,value)
            obj.Lambda = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.PRstruct.Lambda = value;
                obj.PSF4pi.PRstruct1.Lambda = value;
                obj.PSF4pi.PRstruct2.Lambda = value;
            end
        end
        
        function set.RI_imm(obj,value)
            obj.RI_imm = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.PRstruct.RefractiveIndex = value;
                obj.PSF4pi.PRstruct1.RefractiveIndex = value;
                obj.PSF4pi.PRstruct2.RefractiveIndex = value;
            end
        end
        
        function set.RI_med(obj,value)
            obj.RI_med = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.nMed = value;
            end
        end
        
        function set.OTFsigma(obj,value)
            obj.OTFsigma = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.PRstruct.SigmaX = value;
                obj.PSF4pi.psfobj.PRstruct.SigmaY = value;
                obj.PSF4pi.PRstruct1.SigmaX = value;
                obj.PSF4pi.PRstruct1.SigmaY = value;
                obj.PSF4pi.PRstruct2.SigmaX = value;
                obj.PSF4pi.PRstruct2.SigmaY = value;
            end
        end
        
        function set.NA_top(obj,value)
            obj.NA_top = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.PRstruct1.NA = value;
            end
        end
        
        function set.NA_bot(obj,value)
            obj.NA_bot = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.PRstruct2.NA = value;
            end
        end
        
        function set.Iratio(obj,value)
            obj.Iratio = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.Iratio = value;
            end
        end
        
        function set.Modulationdepth(obj,value)
            obj.Modulationdepth = value;
            if ~isempty(obj.PSF4pi)
                obj.PSF4pi.psfobj.ModulationDepth = value;
            end
        end
        
    end
end