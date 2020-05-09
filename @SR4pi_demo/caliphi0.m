function caliphi0(obj,dataname,type)
% caliphi0 - calibrate cavity phase, which is assciated with room temperature and humidity change
%   Data were collected at high laser power at 0.02 s exposure time, one
%   data file was saved for every 2000 frames, one dataset contains 10-90
%   data files. For cavity phase calibration, each data file was segmented
%   into <= 500 frames per segement and one cavity phase was estimated for
%   each segment.
%
% Data format: each data file contains four 3D matrix with dimensions of
%              (x,y,frame) corresponding to data from the four imaging quandrants
%
% Input:
%   dataname - the keyword of the data files to be estimated, e.g. if dataname is 'cell2', then all data files with the keyword 'cell2' under the obj.Datapath will be selected
%   type - 'single': dataset is from single-optical section imaging
%          'multi': dataset is from multi-optical section imaging
%
% Output:
%   ResultPR.Phi0cali - calibration results of cavity phase
%       PhiL - estimated cavity phase for all segments
%       P_phi - estimated slope of the stripe patterns in sr-phi plot for all data files
%       phi_sm - fitted cavity phase with a smooth spline curve and will be used for MLE localization
%
% Output will be stored in obj and saved in a file with the keyword '_phical4pi' under obj.Resultpath\PR fit\
%%
figfolder=[obj.Resultpath,'\PR fit\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end

datafile = dir([obj.Datapath,dataname,'*.mat']);
dirN = length(datafile);

%% prefit phi0
load([obj.Datapath,datafile(1).name])
obj.Maxframe = size(qd1,3);
tN = max([round(obj.Maxframe/500)+1,3]);         % 500 frames per section
ts = linspace(0,obj.Maxframe,tN);
acc_thresh = 700;
plotflag = 0;
Phi0 = zeros(tN-1,1);
PhiL = [];
P_phi = [];
for ii = 1:dirN
    disp(['process file number:',num2str(ii)])
    tic
    % make subregions
    [psf4dt_wT,psf4dt_noT,~,cor,trans_dx,trans_dy,rt] = makerois4pi_sCMOS(obj.Datapath,datafile(ii).name,obj,plotflag);
    % generate initial guess
    ct_thresh = 0.3;
    astfitflag = 0;
    [sr_data,phi_data,cor_data] = genini4pi_loc(psf4dt_wT,psf4dt_noT,cor,ct_thresh,obj,trans_dx,trans_dy,[],astfitflag,plotflag);
    [~,p_phi] = findphi0_loc(sr_data,phi_data,acc_thresh,obj,[],plotflag);
    clear psf4dt_wT psf4dt_noT trans_dx trans_dy
    for nn = 1:tN-1
        mask = cor_data(:,3)>=ts(nn)&cor_data(:,3)<ts(nn+1);
        [phic,p_phi] = findphi0_loc(sr_data(mask),phi_data(mask),acc_thresh,obj,p_phi,plotflag);
        Phi0(nn) = phic;
    end
    toc
    PhiL = cat(1,PhiL,Phi0);
    P_phi = cat(1,P_phi,p_phi);
end

switch type
    case 'single'
        phi_unwrap = unwrap(PhiL);
        v = [1:numel(phi_unwrap)]';
        f_phi0 = fit(v,phi_unwrap,'smoothingspline','SmoothingParam',0.005);
        phi_sm = feval(f_phi0,v);
        figure;plot(v,phi_unwrap,'.')
        hold on
        plot(v,phi_sm,'-')
        xlabel('section number')
        ylabel('\phi_0')
    case 'multi'
        Nt = tN-1;
        Ns = 1+str2double(datafile(dirN).name(end-6:end-4));
        Nc = numel(PhiL)/Nt/Ns;
        nimm = obj.RI_imm;
        nmed = obj.RI_med;
        stepsize = obj.Stepsize; % micron
        lambda = obj.Lambda_bead; % micron
        op = 2*nimm*stepsize*(1-(nmed/nimm)^2);
        phi_step = wrapToPi(op/lambda*2*pi);
        PhiL_v2 = zeros(size(PhiL));
        for s =1:Ns % step index
            ind = [];
            for ii = 1:Nc
                ind = cat(2,ind,1+(s-1)*Nt+Ns*Nt*(ii-1):(s-1)*Nt+Ns*Nt*(ii-1)+Nt);
            end
            PhiL_v2(ind) = PhiL(ind)-(phi_step*(s-1));
        end
        figure;plot(PhiL_v2,'.')
        %
        tol = pi*1; % change this value if necessary
        phi_unwrap_v2 = unwrap(PhiL_v2,tol);
        v = [1:numel(PhiL_v2)]';
        f_phi0 = fit(v,phi_unwrap_v2,'smoothingspline','SmoothingParam',0.0001);
        phi_sm_v2 = feval(f_phi0,v);
        
        if plotflag == 1
            figure;plot(phi_sm_v2,'.')
            hold on;
            plot(unwrap(PhiL_v2,tol),'o')
        end
        %%
        phi_section = [];
        phi_sm = zeros(size(PhiL));
        if plotflag == 1
            figure;
            ha = axes;hold on
        end
        for s =1:Ns % step index
            ind = [];
            for ii = 1:Nc
                ind = cat(2,ind,1+(s-1)*Nt+Ns*Nt*(ii-1):(s-1)*Nt+Ns*Nt*(ii-1)+Nt);
            end
            
            phi_sm(ind) = phi_sm_v2(ind)+(phi_step*(s-1));
            phi_d = unwrap(PhiL(ind))-round(mean(unwrap(PhiL(ind))-phi_sm(ind))./2./pi)*2*pi;
            if plotflag == 1
                plot(ind,phi_sm(ind),'.',ind,phi_d,'o')
            end
            phi_section = cat(2,phi_section,phi_d);
        end
        
        %% refine fit
        phi_step1 = [0,cumsum(diff(mean(phi_section,1)))];
        PhiL_v2 = zeros(size(PhiL));
        for s =1:Ns % step index
            ind = [];
            for ii = 1:Nc
                ind = cat(2,ind,1+(s-1)*Nt+Ns*Nt*(ii-1):(s-1)*Nt+Ns*Nt*(ii-1)+Nt);
            end
            
            PhiL_v2(ind) = PhiL(ind)-phi_step1(s);
        end
        figure;plot(PhiL_v2,'.')
        %
        tol = pi*1; % change this value if necessary
        phi_unwrap_v2 = unwrap(PhiL_v2,tol);
        v = [1:numel(PhiL_v2)]';
        f_phi0 = fit(v,phi_unwrap_v2,'smoothingspline','SmoothingParam',0.0001);
        phi_sm_v2 = feval(f_phi0,v);
        
        figure;plot(phi_sm_v2,'.')
        hold on;
        plot(unwrap(PhiL_v2,tol),'o')
        xlabel('section number')
        ylabel('\phi_0')

        %%
        phi_section = [];
        phi_sm = zeros(size(PhiL));
        figure;
        ha = axes;hold on
        for s =1:Ns % step index
            ind = [];
            for ii = 1:Nc
                ind = cat(2,ind,1+(s-1)*Nt+Ns*Nt*(ii-1):(s-1)*Nt+Ns*Nt*(ii-1)+Nt);
            end
            
            phi_sm(ind) = phi_sm_v2(ind)+phi_step1(s);
            phi_d = unwrap(PhiL(ind))-round(mean(unwrap(PhiL(ind))-phi_sm(ind))./2./pi)*2*pi;
            plot(ind,phi_sm(ind),'.',ind,phi_d,'o')
            phi_section = cat(2,phi_section,phi_d);
        end  
        xlabel('section number')
        ylabel('\phi_0')
end
obj.ResultPR.Phi0cali = struct('PhiL',PhiL,'P_phi',P_phi,'phi_sm',phi_sm);
save([figfolder,datafile(1).name(1:end-4),'_',type,'_phical4pi.mat'],'PhiL','P_phi','phi_sm')

end