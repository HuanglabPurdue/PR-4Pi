function optimparam(obj,type)
% optimparam - 
%   the first time call optimparam(), select the same bead data used for
%   caliparam(). The data is a z-stack of measured 4PiPSFs. Data were collected by
%   imaging a 40 nm fluorescence bead at z positions from -1 to 1 um, with
%   20 nm step size, 0.1 s exposure time and 1 frame per step at low laser power. 
%   Make sure there is no axial drift during data acquisition.
%
% Data format: 4D matrix (x,y,frame,z)
%
% Input:
%   type - 'lambda': optimize emission wavelength (relative to air) for bead data by minimizing axial localization deviation within a defined z range
%          'srcenter': optimize the sr (shape metric) for infocus PSF by minimizing the offset of sr-z curves from data and PSF models
%          'Iratio': optimize intensity ratio between top and bottom emission path by matching the slope of sr-z curves from data and PSF models
%          'evaluate': evaluate the goodness of optimized parameters
%               
% Output:
%   Lambda - emission wave length used for 4PiPSF model, unit: micron 
%   Lambda_bead - emission wave length used for bead data, it is the wavelength relative to air and will be a constant after optimizing the hyper parameters (see optimparam), unit: micron
%   Srcenter - shape metric for infocus astigmatism PSF
%   Iratio - transmission ratio between top and bottom emission path, from 0 to 1 
%   Paramcheck - goodness of optimized parameters: srcenter, lambda, Iratio
%       srcenter_dis: ideally should be 0
%       Iratio_ratio: ideally should be 1, but this value can have larger tolerance
%       lambda_ratio: ideally should be 1
%
% Output will be stored in obj

%%
plotflag = 0;
if isempty(obj.PSFoptidata)
    [datafile, datapath] = uigetfile([obj.Datapath,'*.mat'],'Select file:','MultiSelect','off');
    
    [psf4dt_wT,psf4dt_noT,gainR4dt,cor1,trans_dx,trans_dy,rot] = makerois4pi_sCMOS_bead(datapath,datafile,obj,plotflag);
    obj.PSFoptidata = struct('psf4dt_wT',psf4dt_wT,'psf4dt_noT',psf4dt_noT,'gainR4dt',gainR4dt,'cor',cor1,...
                        'transx',trans_dx,'transy',trans_dy,'rt',rot);
end
Nfit = size(obj.PSFoptidata.transx,1);        % should be 101
Np = 1;
Nz = Nfit/Np;
close all;
switch type
    case 'lambda'
        obj.Lambda = obj.Lambda_bead;
        fitrange = [21:Nz-20];
        [PcudaM] = testparam(obj,Np,plotflag);
        z_pp = PcudaM(:,3).*1e3;
        obj.Lambda = obj.Lambda*(2000/(Nz-1))/mean(diff(z_pp(fitrange)));
        obj.Lambda_bead = obj.Lambda;
    case 'srcenter'
        fitrange = [31:Nz-30];
        src = [-0.2,0.2];
        dis = zeros(1,2);
        for ii = 1:2
            obj.Srcenter = src(ii);
            [~,sr_psf,sr_data] = testparam(obj,Np,plotflag);
            dis(ii) = mean(sr_psf(fitrange)-sr_data(fitrange));
        end
        obj.Srcenter = src(1)-dis(1)*diff(src)/diff(dis);
    case 'Iratio'
        fitrange = [31:Nz-30];
        [~,sr_psf,sr_data] = testparam(obj,Np,plotflag);
        obj.Iratio = obj.Iratio*mean(diff(sr_psf(fitrange)))./mean(diff(sr_data(fitrange)));
    case 'evaluate'
        plotflag = 1;
        [PcudaM,sr_psf,sr_data] = testparam(obj,Np,plotflag);
        fitrange = [31:Nz-30];
        dis = mean(sr_psf(fitrange)-sr_data(fitrange));
        pI = mean(diff(sr_psf(fitrange)))./mean(diff(sr_data(fitrange)));
        z_pp = PcudaM(:,3).*1e3;
        fitrange = [21:Nz-20];
        plambda = (2000/(Nz-1))/mean(diff(z_pp(fitrange)));
        obj.Paramcheck = struct('srcenter_dis',dis,'Iratio_ratio',pI,'lambda_ratio',plambda);
end
end

function [PcudaM,sr_psf,sr_data,z_phifit] = testparam(obj,Np,plotflag)
psf4dt_wT = obj.PSFoptidata.psf4dt_wT;
psf4dt_noT = obj.PSFoptidata.psf4dt_noT;
gainR4dt = obj.PSFoptidata.gainR4dt;
cor1 = obj.PSFoptidata.cor;
trans_dx = obj.PSFoptidata.transx;
trans_dy = obj.PSFoptidata.transy;
rot = obj.PSFoptidata.rt;

ct_thresh = 0;
astfitflag = 0;
[sr_data,phi_data,cor_data,psf4dt_data,x0_data] = genini4pi_loc(psf4dt_wT,psf4dt_noT,cor1,ct_thresh,obj,trans_dx,trans_dy,[],astfitflag,plotflag);
[phic,phic_u] = findphi0_bead_loc(sr_data,phi_data,obj);
zT = obj.PhaseparamM0.zT;
z_phifit = (unwrap(phi_data)-phic_u).*zT./2./pi;
%%
[PcudaM,tmp_data,tmp_psf] = locbead(obj,phic,cor1,gainR4dt,x0_data,psf4dt_noT,trans_dx,trans_dy,rot,Np);
[sr_psf,sr_data] = getmetric(tmp_psf,tmp_data);
%% plot metric
if plotflag == 1
    h = figure;
    h.Position = [100,200,500,400];
    ha = axes;hold on;
    plot(PcudaM(:,3),sr_data,'.',PcudaM(:,3),sr_psf,'.')
    ha.XLabel.String = 'z postion (um)';
    ha.YLabel.String = 'shape metric';
    axis tight
    legend('data','model')
end
%% bias plot
pixelsize = obj.Pixelsize;
Nfit = numel(phi_data);
Nz = Nfit/Np;
z_stage = linspace(-1000,1000,Nz);
z_true = reshape(repmat(z_stage,Np,1),Nfit,1);
z_ct = z_phifit-median(z_phifit-z_true);
z_ct = z_ct + round((z_true-z_ct)./zT).*zT;
x_ct = (x0_data(:,1)-median(x0_data(:,1))).*pixelsize.*1e3;
y_ct = (x0_data(:,2)-median(x0_data(:,2))).*pixelsize.*1e3;

z_pp = PcudaM(:,3).*1e3-median(PcudaM(:,3).*1e3-z_true);
x_pp = (PcudaM(:,1)-median(PcudaM(:,1))).*pixelsize.*1e3;
y_pp = (PcudaM(:,2)-median(PcudaM(:,2))).*pixelsize.*1e3;

if plotflag == 1
    lw = 1;
    poslabel = {'x','y','z'};
    xLim = [-1000,1000];
    yLim = [-50,50];
    h2 = figure;
    h2.Position = [150,150,470,630];
    for ii = 1:3
        ha3(ii) = subplot(3,1,ii);
        if ii==3
            pos_true = z_true;
        else
            pos_true = 0;
        end
        eval(['pos=',poslabel{ii},'_ct;']);
        ct_bias = mean(reshape(pos-pos_true,Np,Nz),1);
        eval(['pos=',poslabel{ii},'_pp;']);
        pp_bias = mean(reshape(pos-pos_true,Np,Nz),1);
        plot(z_stage,ct_bias,'-r^',z_stage,pp_bias,'-bo','linewidth',1.2)
        hold on;
        line([xLim],[20,20],'color',[0,0,0],'linestyle','--','linewidth',lw);
        line([xLim],[-20,-20],'color',[0,0,0],'linestyle','--','linewidth',lw);
        ha3(ii).XLim = xLim;
        ha3(ii).YLim = yLim;
        ha3(ii).FontSize = 12;
        ha3(ii).XLabel.String = ['stage position (nm)'];
        ha3(ii).YLabel.String = [poslabel{ii},' deviation (nm)'];
        
    end
end
end

function [sr1,sr2] = getmetric(tmp_psf,tmp_data)
PSFsigma = 1.15;
iterations = 100;
fittype = 4;
datasum = squeeze(sum(tmp_psf,4));
datasum(datasum<=0) = 1e-6;
[P, CRLB, LL] = GPUgaussMLEv2(single(datasum),PSFsigma,iterations,fittype);

sx1 = P(:,5);
sy1 = P(:,6);
sr1 = sx1.^2 - sy1.^2;

datasum = squeeze(sum(tmp_data,4));
datasum(datasum<=0) = 1e-6;
[P, CRLB, LL] = GPUgaussMLEv2(single(datasum),PSFsigma,iterations,fittype);

sx2 = P(:,5);
sy2 = P(:,6);
sr2 = sx2.^2 - sy2.^2;

end

function  [PcudaM,tmp_data,tmp_psf] = locbead(obj,phic,cor1,gainR4dt,x0_data,psf4dt_noT,trans_dx,trans_dy,rot,Np)
plotflag = 0;
psfsize = 256;
boxsize = 25;
bin = 4;
Nzs = 301;
[samplepsf,startx,starty,startz,sizex,dz,dx] = gensamplepsf4pi(obj,psfsize,boxsize,phic,bin,Nzs);

samplepsf_4d = [];
samplepsf_cuda = [];

for ii = 1:4
    samplepsf_4d = cat(4,samplepsf_4d,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
    samplepsf_cuda = cat(3,samplepsf_cuda,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
end
samplepsf_cuda = single(samplepsf_cuda);

[st] = genpsfstruct(samplepsf_4d(:,:,:,1),dx,dz,'matrix');
for ii = 2:4
    [sti] = genpsfstruct(samplepsf_4d(:,:,:,ii),dx,dz,'matrix');
    st = catstruct(st,sti,3);
end
sobj = struct('samplepsf_cuda',samplepsf_cuda,'dx',dx,'dz',dz,'startx',startx,'starty',starty,'startz',startz);
obj.PSFlib = sobj;
%%
data0 = psf4dt_noT;
data0(data0<=0) = 1e-6;
[z_guess,mcc_max] = genini4pi_z_cuda(data0,sobj);
if plotflag == 1
    figure;plot(z_guess.*1e3,'.')
    figure;plot(mcc_max)
end
Nfit = size(trans_dx,1);
Nz = Nfit/Np;
z_stage = linspace(-1000,1000,Nz);
z_true = reshape(repmat(z_stage,Np,1),Nfit,1);
z_guess = z_true./1e3;

%%
data = psf4dt_noT;
data(data<=0) = 1e-6;
Nfit = size(data,3);
boxsize_data = obj.Boxsize;
data_cuda = single(data(:)); %psfsize x Nfit x quadrantN
coords_cuda = single(cor1(:));
gainRsub = single(gainR4dt(:));

xtmp = x0_data(:,1)+boxsize_data/2;
ytmp = x0_data(:,2)+boxsize_data/2;
x0 = cat(2,xtmp,ytmp,z_guess,x0_data(:,4:11));
x0i = single(reshape(x0',Nfit*11,1));
lambda0 = 1;% damping factor
tic
iterateN = 100;
[P_cuda,CG,crlb,err,PSF_cuda] = cuda4pi_loc_spline(data_cuda,coords_cuda,single(st.F),...
                                dx,dz,startx,starty,startz,iterateN,Nfit,lambda0,x0i,gainRsub,...
                                single(st.Fx),single(st.Fy),single(st.Fz),...
                                single(st.Fxy),single(st.Fxz),single(st.Fyz),single(st.Fxyz),...
                                single(trans_dx(:)),single(trans_dy(:)),single(rot));


toc
clear cuda4pi_loc_spline

% check result
tmp_psf = reshape(PSF_cuda,[boxsize_data,boxsize_data,Nfit,4]);
tmp_data =reshape(data_cuda,[boxsize_data,boxsize_data,Nfit,4]);
psf_cuda = [];
Data_cuda = [];
for ii = 1:4
psf_cuda = cat(2,psf_cuda,tmp_psf(:,:,:,ii));
Data_cuda = cat(2,Data_cuda,tmp_data(:,:,:,ii));
end
ov = joinchannels('RGB',Data_cuda,psf_cuda);
%dipshow(ov)
%dipshow(cat(1,Data_cuda,psf_cuda))
crlbM = reshape(crlb,11,Nfit)';
stdM = sqrt(crlbM);
errM = reshape(err,2,Nfit)';
cgM = reshape(CG,11,Nfit)';
PcudaM = reshape(P_cuda,11,Nfit)';

end