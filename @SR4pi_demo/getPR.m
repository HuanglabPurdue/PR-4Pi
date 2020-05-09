function getPR(obj,type)
% getPR - generate phase retrieved pupil functions
%   select a z-stack of measured PSFs from either bottom or top emission path.
%   Data were collected by imaging a 40 nm fluorescence bead at z positions from -1 to 1 um, with
%   100 nm step size, 1 s exposure time and 1 frame per step at low laser power. 
%   Make sure there is no axial drift during data acquisition.
%
% Data format: 4D matrix, (x,y,frame,z)
%
% Input:
%   type - 'bottom': phase retrieval for bottom emission path
%          'top': phase retrieval for top emission path
%
% Output:
%   PSFtop - object of OptimPR_Ast class for phase retrieval of top emission path
%   PSFbot - object of OptimPR_Ast class for phase retrieval of bottom emission path
%
% Output will be stored in obj and saved in file with key word '_result_optimAst_botobj' or 
% '_result_optimAst_topobj' under obj.Resultpath\PSF measure\PR_result\
%%
[FileName, FileDir] = uigetfile([obj.Datapath,'*.mat'],'Select file:','MultiSelect','off');

oprobj = OptimPR_Ast();
oprobj.PRobj.CCDoffset = obj.Constantoffset;
oprobj.PRobj.Gain = obj.Constantgain;
switch type
    case 'bottom'
        oprobj.PRobj.PRstruct.NA = obj.NA_bot;
        oprobj.PRobj.Zstart = -1;           %micron
        oprobj.PRobj.Zend = 1;              %micron
        oprobj.PRobj.Zstep = 0.1;           %micron
    case 'top'
        oprobj.PRobj.PRstruct.NA = obj.NA_top;
        oprobj.PRobj.Zstart = 1;            %micron
        oprobj.PRobj.Zend = -1;             %micron
        oprobj.PRobj.Zstep = -0.1;          %micron
end
oprobj.PRobj.PRstruct.Lambda = obj.Lambda;
oprobj.PRobj.PRstruct.RefractiveIndex = obj.RI_imm;
oprobj.PRobj.nMed = obj.RI_med;
oprobj.PRobj.Pixelsize = obj.Pixelsize;
oprobj.PRobj.PRstruct.SigmaX = obj.OTFsigma;
oprobj.PRobj.PRstruct.SigmaY = obj.OTFsigma;
oprobj.FileDir = FileDir;
oprobj.FileName = FileName;
% the followings parameters are fixed 
oprobj.PRobj.PSFsize = 128;
oprobj.PRobj.SubroiSize = 34;               % we cut the initial image. Different sizes might give better or worse results.
oprobj.PRobj.OTFratioSize = 60;
oprobj.PRobj.ZernikeorderN = 7;
oprobj.PRobj.Zindstart = 1;                 %index
oprobj.PRobj.Zindend = 21;
oprobj.PRobj.Zindstep = 4;
oprobj.PRobj.IterationNum = 25;
oprobj.PRobj.IterationNumK = 5;
oprobj.PRobj.mvType = 'mvstage';
oprobj.PRobj.PSFtype = 'IMM';
oprobj.PRobj.Stagepos = 0; %micron
oprobj.PRobj.Ztype = 'uniform';
oprobj.FitZrange = [-0.8,0.8];
oprobj.IterationMonte = 50;
oprobj.PRobj.Enableunwrap = 0;

%% generate initial PR result
oprobj.prepdata('EMCCD');
oprobj.initialPR();
oprobj.PRobj.genPRfigs('zernike');
oprobj.PRobj.genPRfigs('pupil');
oprobj.PRobj.genPRfigs('PSF');
%oprobj.PRobj.calcrlb();
%% save OptimPR_Ast object
resdir=fullfile(obj.Resultpath,'\PSF measure\PR_result\');
if ~exist(resdir,'dir')
    mkdir(resdir)
    
end
oprobj.SaveDir=resdir;
switch type
    case 'bottom'
        obj.PSFbot = oprobj;
        oprobj.SaveName='_result_optimAst_botobj';
    case 'top'
        obj.PSFtop = oprobj;
        oprobj.SaveName='_result_optimAst_topobj';
end
oprobj.saveObj();


end