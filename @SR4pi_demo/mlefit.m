function mlefit(obj,dataname)
% mlefit - localization using maximum likelihood estimation (MLE) based on PR-4PiPSF model
%   Data were collected at high laser power at 0.02 s exposure time, one
%   data file was saved for every 2000 frames, one dataset contains 10-90
%   data files. For localization, each data file was segmented into <= 500 frames  
%   per segement and one PR-4PiPSF library (PSFlib) was generated for each segment. 
%   The PSFlib covers an axial range of -1.5 to 1.5 µm and a lateral 
%   extend of 25 × 25 camera pixels, and each voxel of the PSFlib measures 
%   0.25 pixel × 0.25 pixel × 10 nm. Derivatives and PSF model used in MLE
%   were generated through spline interpolation of the PSFlib and the cost
%   function was minimized through the Levenberg-Marquardt algorithm. 
%
%   mlefit() includes the following major procedures:
%       initial estimation of x,y,I,bg
%       initial estimation of z
%       first MLE localization (on GPU)
%       ghost reduction
%       section MLE localization (on GPU)
%
% Data format: each data file contains four 3D matrix with dimensions of
%              (x,y,frame) corresponding to data from the four imaging quandrants
%
% Input:
%   dataname - the keyword of the data files to be estimated, e.g. if dataname is 'cell2', 
%              then all data files with the keyword 'cell2' under the obj.Datapath will be selected
%
% Output:
%   ResultPR.MLE - MLE localization results from PR-4Pi algorithm
%       res - structure of localization results
%           x - x positions, unit: nm
%           y - y positions, unit: nm
%           z - z positions, unit: nm
%           zast - z positions from astigmatism fitting with a Gaussian PSF model, unit: nm
%           I - photon per objective
%           bg - background per pixel
%           LL - loglikelihood ratio
%           STD - theoretical localization precisions from CRLB for x (um),y (um),z (um),I,bg
%           phi - interference phase
%           sx - width of the astigmatism PSFs along the x-axis
%           sy - width of the astigmatism PSFs along the y-axis
%           cg - convergence of MLE
%           err - SSE and LLR of each localization
%           mcc - maximum cross correlation value from initial estimation of z
%           f - data file number
%           t - frame number relative to each data file
%           stepN - optical section number
%           cycN - cycle number of imaging through one cycle of all optical sections, it is the same as f for single-optical section imaging
%       Phifit - estimated modulation period (nm) and the cavity phase for each data segment from localization results
%       gper - ratio of ghost pairs after and before ghost reduction for each data segment
%
% Output will be stored in obj and saved in a file with the keyword '_loc4pi_constI' under obj.Resultpath\PRfit\

%%
figfolder=[obj.Resultpath,'\PR fit\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end

datafile = dir([obj.Datapath,dataname,'*.mat']);
dirN = length(datafile);

%%
tN = max([round(obj.Maxframe/500)+1,3]);         % 500 frames per section
ts = linspace(0,obj.Maxframe,tN);
plotflag = 0;
Res = struct('transx',[],'transy',[],'sx',[],'sy',[],'x0',[],'cor',[],'phi',[],'zast',[],'z_guess',[],'mcc',[],'P',[],'err',[],...
    'cg',[],'crlb',[],'STD',[]);
filenum = [];
cycnum = [];
stepnum = [];
gper = [];
Phifit = [];
P_phi = obj.ResultPR.Phi0cali.P_phi;
pixelsize = obj.Pixelsize;
p_phi = median(P_phi);
phi_sm = obj.ResultPR.Phi0cali.phi_sm;
mccflag = 0;
obj.Lambda = obj.Lambda_bead*obj.RI_imm/obj.RI_med;
%obj.Lambda = obj.Lambda_bead;
for ii = 1:dirN
    %p_phi = P_phi(ii);
    disp(['process file number:',num2str(ii)])
    tic
    % make subregions
    [psf4dt_wT,psf4dt_noT,gainR4dt,cor,trans_dx,trans_dy,rt] = makerois4pi_sCMOS(obj.Datapath,datafile(ii).name,obj,plotflag);
    % generate initial guess
    ct_thresh = 0.3;
    astfitflag = 1;
    [sr_data,phi_data,cor_data,psf4dt_data,x0_data,~,z_ast,sx_data,sy_data,maskall] = genini4pi_loc(psf4dt_wT,psf4dt_noT,cor,ct_thresh,obj,trans_dx,trans_dy,p_phi,astfitflag,plotflag);
    
    gainR4dt_data = gainR4dt(:,:,maskall,:);
    transx_data = trans_dx(maskall,:);
    transy_data = trans_dy(maskall,:);
    
    res0 = struct('transx',transx_data,'transy',transy_data,'sx',sx_data,'sy',sy_data,...
        'x0',x0_data,'cor',cor_data,'phi',phi_data,'zast',z_ast);
    
    clear psf4dt_wT psf4dt_noT gainR4dt cor trans_dx trans_dy transx_data transy_data sx_data sy_data x0_data cor_data
    time1 = toc;
    disp(['Elapsed time for subregion selection: ',num2str(time1),' s'])
    gper1 = [];
    phi1 = [];
    resf = struct('transx',[],'transy',[],'sx',[],'sy',[],'x0',[],'cor',[],'phi',[],'zast',[],'z_guess',[],'mcc',[],'P',[],'err',[],...
    'cg',[],'crlb',[],'STD',[]);

    for nn = 1:tN-1
        
        mask = res0.cor(:,3)>=ts(nn)&res0.cor(:,3)<ts(nn+1);
        gainR4dt_sub = gainR4dt_data(:,:,mask,:);
        psf4dt_sub = psf4dt_data(:,:,mask,:);
        res1 = applymask(res0,mask);
        %% generate sample psf
        tic
        phic = phi_sm((ii-1)*(tN-1)+nn);
        bin = 4;
        Nzs = 301;
        psfsize = 256;
        boxsize = 25;
        
        [samplepsf,startx,starty,startz,sizex,dz,dx] = gensamplepsf4pi(obj,psfsize,boxsize,phic,bin,Nzs);
        samplepsf_4d = [];
        samplepsf_cuda = [];
        for ss = 1:4
            samplepsf_4d = cat(4,samplepsf_4d,permute(reshape(samplepsf{ss},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
            samplepsf_cuda = cat(3,samplepsf_cuda,permute(reshape(samplepsf{ss},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
        end
        samplepsf_cuda = single(samplepsf_cuda);
        [st] = genpsfstruct(samplepsf_4d(:,:,:,1),dx,dz,'matrix');
        for ss = 2:4
            [sti] = genpsfstruct(samplepsf_4d(:,:,:,ss),dx,dz,'matrix');
            st = catstruct(st,sti,3);
        end
        sobj = struct('samplepsf_cuda',samplepsf_cuda,'st',st,'dx',dx,'dz',dz,'startx',startx,'starty',starty,'startz',startz);
        obj.PSFlib = sobj;
        clear st sti samplepsf_cuda
        time1 = toc;
        disp(['Elapsed time for library 4Pi-PSF generation: ',num2str(time1),' s'])
        %% initial guess in z
        data0 = psf4dt_sub;
        data0(data0<=0) = 1e-6;
        [z_guess,mcc_max] = genini4pi_z_cuda(data0,sobj);
        res1.z_guess = z_guess;
        res1.mcc = mcc_max;
        if plotflag == 1
            figure;hist(z_guess,100)
            figure;hist(mcc_max,100)
        end
        if mccflag == 1
            maskc = mcc_max>0.6;
        else
            maskc = mcc_max>0;
        end
        gainR4dt_mcc = gainR4dt_sub(:,:,maskc,:);
        psf4dt_mcc = psf4dt_sub(:,:,maskc,:);
        res2 = applymask(res1,maskc);
        Nfit = sum(maskc);
        %% localization spline interpolation
        data = psf4dt_mcc;
        data(data<=0) = 1e-6;
        [PcudaM,errM] = loc4pi_cuda_sCMOS_constI(data,gainR4dt_mcc,obj,res2,rt,plotflag);
        res2.P = PcudaM;
        res2.err = errM;
        
        if plotflag ==1
            figure;plot(res2.sx.^2-res2.sy.^2,PcudaM(:,3),'.');
        end
        %% remove ghost
        [~,zT_fit,phi0_fit] = rmoutlier_phi(res2.phi,PcudaM(:,3),phic,plotflag);
        phi_rm = [zT_fit,phi0_fit];
        gzone.xlim = [0,50];
        gzone.zlim = [-25,25]+zT_fit;
        gzone.zT = zT_fit;
        tic
        [~,gpair,z_guess1] = getgpair(res2,gzone,pixelsize,plotflag);
        time1 = toc;
        disp(['Elapsed time for ghost reduction: ',num2str(time1),' s'])

        %% refine z position
        res3 = res2;
        res3.x0 = PcudaM;
        res3.x0(:,1) = PcudaM(:,1)-obj.Boxsize/2;
        res3.x0(:,2) = PcudaM(:,2)-obj.Boxsize/2;
        res3.z_guess = z_guess1./1e3;
        
        [PcudaM,errM,cgM,crlbM,STDM] = loc4pi_cuda_sCMOS_constI(data,gainR4dt_mcc,obj,res3,rt,plotflag);
        res3.P = PcudaM;
        res3.err = errM;
        res3.cg = cgM;
        res3.crlb = crlbM;
        res3.STD = STDM;
        tic
        [~,gpair1] = getgpair(res3,gzone,pixelsize,plotflag);
        time1 = toc;
        disp(['Elapsed time for ghost reduction: ',num2str(time1),' s'])
        
        if plotflag == 1
            sr_fit = res3.sx.^2-res3.sy.^2;
            figure;plot(sr_fit,res2.P(:,3).*1e3,'.')
            figure;plot(sr_fit,res3.P(:,3).*1e3,'.')
        end
        
        if plotflag == 1
            x_fit = (res2.P(:,1)+res2.cor(:,1)).*pixelsize.*1e3;
            y_fit = (res2.P(:,2)+res2.cor(:,2)).*pixelsize.*1e3;
            z_fit = res2.P(:,3).*1e3;
            
            x_fit1 = (res3.P(:,1)+res3.cor(:,1)).*pixelsize.*1e3;
            y_fit1 = (res3.P(:,2)+res3.cor(:,2)).*pixelsize.*1e3;
            z_fit1 = res3.P(:,3).*1e3;
            
            figure;plot(y_fit,z_fit,'.');figure;plot(y_fit1,z_fit1,'.')
            figure;plot(x_fit,z_fit,'.');figure;plot(x_fit1,z_fit1,'.')
        end
        
        resf = catstruct(resf,res3,1);
        phi1 = cat(1,phi1,phi_rm);
        
        gper1 = cat(1,gper1,length(gpair1)/length(gpair));
        disp(['number of localization: ',num2str(Nfit)])
    end
    % concatenate result
    Res = catstruct(Res,resf,1);
    Nffit = numel(resf.P(:,1));
    filenum = cat(1,filenum,ii.*ones(Nffit,1));
    stepnum = cat(1,stepnum,1+ones(Nffit,1).*str2double(datafile(ii).name(end-6:end-4)));
    cycnum = cat(1,cycnum,1+ones(Nffit,1).*str2double(datafile(ii).name(end-10:end-8)));
    gper = cat(1,gper,gper1);
    Phifit = cat(1,Phifit,phi1);
end
Res.f = filenum;
Res.stepN = stepnum;
Res.cycN = cycnum;
Res.I = Res.P(:,4);
Res.bg = Res.P(:,5);
Res.LL = Res.err(:,2);
Res.t = Res.cor(:,3);
x_fit = (Res.P(:,1)+Res.cor(:,1)-obj.Boxsize/2-0.5).*pixelsize.*1e3;
y_fit = (Res.P(:,2)+Res.cor(:,2)-obj.Boxsize/2-0.5).*pixelsize.*1e3;
z_fit = Res.P(:,3).*1e3;
Res.x = x_fit;
Res.y = y_fit;
Res.z = z_fit;
Res.zast = Res.zast.*1e3;        
res = rmfield(Res,{'P','cor'});
%res.z = res.z*obj.RI_imm/obj.RI_med;
close(gcf)
save([figfolder,datafile(1).name(1:end-4),'_loc4pi_constI.mat'],'res','Phifit','gper')
obj.ResultPR.MLE = struct('res',res,'Phifit',Phifit,'gper',gper);
obj.Savename = datafile(1).name(1:end-4);
end