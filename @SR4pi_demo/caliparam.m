function caliparam(obj)
% caliparam - calibrate parameters of incoherent (astigmatism) and interferometric PSFs.
%   select a z-stack of measured 4PiPSFs. Data were collected by
%   imaging a 40 nm fluorescence bead at z positions from -1 to 1 um, with
%   20 nm step size, 0.1 s exposure time and 1 frame per step at low laser power. 
%   Make sure there is no axial drift during data acquisition.
%   
% Data format: 4D matrix (x,y,frame,z)
%               
% Output:
%   PSFcalidata - 4PiPSF data used for calibration, including a z-stack PSFs for each quadrant, see caliparam
%   Astparamc - calibration parameters for Gaussian fitting of astigmatism PSFs, used in contrast method
%   Astparamp - calibration parameters for Gaussian fitting of astigmatism PSFs, used in PR-4Pi method
%   PhaseparamM0 - calibration parameters of the central moment modulation of the 4PiPSFs
%   PhaseparamM3 - calibration parameters of the third moment modulation of the 4PiPSFs
%
% Output will be stored in obj
%%
figfolder=[obj.Resultpath,'\PSF measure\'];
if ~exist(figfolder,'dir')
mkdir(figfolder)
end

[FileName, FileDir] = uigetfile([obj.Datapath,'*.mat'],'Select file:','MultiSelect','off');

load(fullfile(FileDir,FileName))
%%
plotflag = 0;
for ii = 1:4
    eval(['qd = qd',num2str(ii),';'])
    aveqd = squeeze(mean((qd-obj.Constantoffset)./obj.Constantgain,3));
    eval(['aveqd',num2str(ii),' = aveqd;'])
end
clear('qd1','qd2','qd3','qd4')

[q1, q2, q3, q4] = align(aveqd1,aveqd2,aveqd3,aveqd4,obj.Affs);
im1 = cat(2,q3,q4);
im2 = cat(2,q1,q2);
colorim=joinchannels('RGB',im1,im2);
if plotflag == 1
    h = dipshow(colorim);
    diptruesize(h,400)
end
% find local centers
sumim1 = q1+q4+q2+q3;
subsz = 21;
[~,centers] = pickbead(sumim1,subsz,1,'stack',[]);
qd1 = pickbead(q1,subsz,1,'stack',centers);
qd2 = pickbead(q2,subsz,1,'stack',centers);
qd3 = pickbead(q3,subsz,1,'stack',centers);
qd4 = pickbead(q4,subsz,1,'stack',centers);
Nz = size(qd1,3);
%close(gcf)
if plotflag == 1
    h = dipshow(cat(2,qd1,qd2,qd3,qd4));
    diptruesize(h,1000)
end
save([figfolder, '\', FileName(1:end-4), '_psf_' datestr(now,'mm-dd-yy')],'qd1','qd2','qd3','qd4');
obj.PSFcalidata = struct('qd1',qd1,'qd2',qd2,'qd3',qd3,'qd4',qd4); 
%% fit sigmax sigmay
cutim = qd1+qd2+qd3+qd4;
P = GPUgaussMLEv2(single(cutim),1.4,200,4);
sx = P(:,5);
sy = P(:,6);
sr = sx.^2 - sy.^2;

%% generate ast calibration parameter for contrast fitting
fitrange = [10:Nz-10];
stepsize = 20;%nm
[estx,esty,sigmax,sigmay,zstep] = genastparam(P,stepsize,fitrange,figfolder,FileName);
obj.Astparamc = struct('estx',estx,'esty',esty,'zstep',zstep,'sx',sigmax,'sy',sigmay);

%% generate phase calibration parameter
subims=cat(4,qd1,qd2,qd3,qd4);
[phi_s,phi_p,ang_coeff,~,~,phi] = genphaseparam(P,subims,plotflag);
zT = 2*pi/ang_coeff(1);
obj.PhaseparamM0 = struct('phis',phi_s,'phip',phi_p,'ang_coeff',ang_coeff,'zT',zT,'phi',phi);

%% third central moment
subims=cat(4,qd1,qd2,qd3,qd4);
M = 3;
[phi_s,phi_p,ang_coeff,~,~,phi] = genphaseparamM(P,subims,M,plotflag);
zT = 2*pi/ang_coeff(1);
obj.PhaseparamM3 = struct('phis',phi_s,'phip',phi_p,'ang_coeff',ang_coeff,'zT',zT,'phi',phi);

%% find p_phi 
phia = unwrap(obj.PhaseparamM0.phi);
fitrange = [31:Nz-30];
p_phi = polyfit(sr(fitrange),phia(fitrange),1);
f1 = polyval(p_phi,sr);
if plotflag == 1
    figure;
    plot(sr(fitrange),phia(fitrange),'.')
    hold on
    plot(sr,f1)
    xlabel('sigma metric')
    ylabel('interference phase')
end
%% generate ast calibration parameter for pupil fitting
fitrange = [5:Nz-1];
zpos0 = linspace(-1,1,Nz);                      % micron
Ax = -0.1;               
Bx = -0.1;
Ay = 0.1;
By = 0.1;
gamma = 0.45 ;                                  % separation of x/y focal planes
d = 0.5;
PSFsigmax = 1.2;

options = optimset('MaxFunEvals',10000,'MaxIter',10000);
startpoint = [Ax,Bx,Ay,By,gamma,d,PSFsigmax];
estx = fminsearch(@(x) astfuneval(x,zpos0(fitrange)',sx(fitrange),sy(fitrange),2),startpoint,options);

startpoint = [Ax,Bx,Ay,By,-gamma,d,PSFsigmax];
esty = fminsearch(@(x) astfuneval(x,zpos0(fitrange)',sx(fitrange),sy(fitrange),3),startpoint,options);

[~,sx_fit,~]=astfuneval(estx,zpos0');
[~,~,sy_fit]=astfuneval(esty,zpos0');
if plotflag == 1
    figure('position',[200,300,500,400],'color',[1,1,1])
    plot(zpos0,sx,'ro',zpos0,sy,'bo','markersize',5)
    hold on
    plot(zpos0,sx_fit,'r-',zpos0,sy_fit,'b-','linewidth',2)
    legend('Obs. \sigma_x','Obs. \sigma_y','Fit. \sigma_x','Fit. \sigma_y','Location','North');
    xlabel('stage position (um)')
    ylabel('sigma')
end
%% find z estimator
fitrange = [15:Nz-15];
p = polyfit(sr(fitrange)',zpos0(fitrange),1);
f = polyval(p,sr);
if plotflag == 1
    figure;plot(sr,zpos0,'o',sr,f,'-')
    xlabel('sigma metric')
    ylabel('stage position (um)')
end
%% astigmatism fitting with Gaussian PSF
astparam.estx = estx;
astparam.esty = esty;
parfor ii=1:length(sx)
    [zf(ii), zerr(ii)] = iniz_astfit(double(sx(ii)),double(sy(ii)),astparam,0.1);
end

if plotflag == 1
    figure;
    plot(zpos0,zf,'.')
    xlabel('stage position (um)')
    ylabel('found z position from ast. fitting (um)')
end
obj.Astparamp = struct('estx',estx,'esty',esty,'sx',sx,'sy',sy,'sr',sr,'p_phi',p_phi,'p',p);
end


