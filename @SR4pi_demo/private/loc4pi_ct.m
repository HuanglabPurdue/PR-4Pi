% for v14 iPALMast_analysisv14scmos.m

function [Res] = loc4pi_ct(datapath,datafiles,resfolder,obj)
scmos_cali_file = obj.Gainpath;
centers = obj.Quadcenter;
affS = obj.Affs;
thresh = obj.Peakthresh;
llthresh = 800;
Ithresh = 100;
%% parameters
warning('off','all')
fn = length(datafiles);
tmpld = load(scmos_cali_file);
offsetim = tmpld.ccdoffset; % for scmos, yes for now
varim = tmpld.ccdvar+0.1;
gainim = tmpld.gain;
caliims = cat(3,offsetim,varim,gainim);
Res = struct('x',[],'y',[],'zast',[],'phi',[],'t',[],'zast_err',[],'LL',[],'I',[],'stepN',[],'cycN',[],'bg',[],'sx',[],'sy',[],'ct',[],...
    'f',[],'rms',[],'rmp',[]);
%% load one quadrant of the stack
for ff = 1:fn
    %close all
    disp(['process file number:',num2str(ff)])
    tic
    if iscell(datafiles)
        filename = datafiles{ff};
    else
        filename = datafiles(ff).name;
    end
    qds = load(fullfile(datapath,filename));
    flipsigns = [0 0 1 1];
    [calicrops] = iPALMscmos_makeqds(caliims,centers,flipsigns);

    tmpim = calicrops(:,:,3,:);
    tmpim(tmpim<1.3|tmpim>3.5) = mean(mean(mean(calicrops(:,:,3,:))));
    calicrops(:,:,3,:) = tmpim;
    offsetim = squeeze(calicrops(:,:,1,:));
    varim = squeeze(calicrops(:,:,2,:));
    gainim = squeeze(calicrops(:,:,3,:));
    dim = [1 1 size(qds.qd1,3) 1];
    for ii =1:4
            ims = (qds.(['qd',num2str(ii)])-repmat(offsetim(:,:,ii),dim))./repmat(gainim(:,:,ii),dim);
            eval(['qd',num2str(ii),'=ims;']);
    end
    
%     img1 = cat(2,qd1,qd3);
%     img2 = cat(2,qd4,qd2);
%     img_4qds = cat(1,img1,img2);
%     dipshow(img_4qds);
    %clear mex
    %% rotate and align
    [q1, q2, q3, q4] = align(qd1,qd2,qd3,qd4,affS);
    [v1, v2, v3, v4] = align(varim(:,:,1),varim(:,:,2),varim(:,:,3),varim(:,:,4),affS);
    [g1, g2, g3, g4] = align(gainim(:,:,1),gainim(:,:,2),gainim(:,:,3),gainim(:,:,4),affS);
    maskg = g1<=1|g2<=1|g3<=1|g4<=1;
    v1(maskg) = 1e7;
    v2(maskg) = 1e7;
    v3(maskg) = 1e7;
    v4(maskg) = 1e7;
    
    g1(maskg) = 1e-7;
    g2(maskg) = 1e-7;
    g3(maskg) = 1e-7;
    g4(maskg) = 1e-7;
    
    %ov1 = joinchannels('RGB',q1,q3);
    %ov2 = joinchannels('RGB',q2,q4);
    %% find local centers
    sumim1 = q1+q2+q3+q4;
    sumv = v1+v2+v3+v4;
    sumvg = v1./g1./g1+v2./g2./g2+v3./g3./g3+v4./g4./g4;
    sumim1(sumim1<=1e-6) = 1e-6;
    if ff == 1
        h = dipshow(mean(sumim1,3));
        diptruesize(h,400)
        colormap(grey)
        print(gcf,'-dpng','-r300',[resfolder,filename(1:end-4),'_sumprojection_image'])
    end
    
    subsz = obj.Boxsize-1;   % sub region size
    pick_flag = 1; % use range
    [sub_regions, tlz, locmaxc] = iPALMast_sumim_seg(sumim1,thresh,subsz,sumv,sumvg,pick_flag);
    %% fit sigmax sigmay
    disp(['A total of ' num2str(size(sub_regions,3)) ' subregions were detected. Start sCMOS_sigmaxy fitting']);
    PSFsigma = 1.15;
    Iterations = 100;
    FitType = 4;
    
    N = size(sub_regions,3);
    maxN = 10000;
    v0 = [1:maxN:N];
    Ns = numel(v0);
    if v0(Ns)<N
        v0(Ns+1) = N;
    end
    v0(end) = v0(end)+1;
    
    P = [];
    CRLB = [];
    LL = [];
    
    for ii = 1:numel(v0)-1
        vi = v0(ii):v0(ii+1)-1;
        [P1, CRLB1, LL1] = GPUgaussMLEv2_sCMOS(single(sub_regions(:,:,vi)),single(tlz(vi,1:2)'),single(sumvg),PSFsigma,Iterations,FitType);
        
        P = cat(1,P,P1);
        CRLB = cat(1,CRLB,CRLB1);
        LL = cat(1,LL,LL1);
    end
    
    res = struct('x',P(:,2),'y',P(:,1),'I',P(:,3),'bg',P(:,4),'LL',-2*LL,'sx',P(:,5),'sy',P(:,6),...
        'cor',tlz);
    %% filter
    maskxy = res.x<subsz-1 & res.x>1 & res.y<subsz-1 & res.y>1;
    masks = res.sx<subsz/2 & res.sy<subsz/2;
    maskll = res.LL<llthresh; 
    maskothers = res.I>Ithresh;
    
    mask = maskxy&masks&maskll&maskothers;
    res1 = applymask(res,mask);
    locmaxc_f = locmaxc(mask,:);
    res1.stepN = ones(numel(res1.x),1).*str2double(filename(end-6:end-4))+1;
    res1.cycN = ones(numel(res1.x),1).*str2double(filename(end-10:end-8))+1;
    
    %% z ast initial guess
    zf=[];
    zerr=[];
    
    astparam = obj.Astparamc;
    astparam.model = [];
    sx = res1.sx;
    sy = res1.sy;
    parfor ii = 1:length(res1.x)
        [zf(ii), zerr(ii)] = iPALMast_iniz_astfit(sx(ii),sy(ii),astparam);
    end
    res1.zast = zf';
    res1.zast_err = zerr';
    
    %% determine the phase
    
    if isempty(locmaxc_f)
        continue
    end
    [z_ang,rms,rmp,ang_ctr] = iPALMast_RM_zang(obj,qd1,qd2,qd3,qd4,subsz,varim,gainim,res1);
    
    %% collect data
    pixelsize = obj.Pixelsize*1e3;
    res1.x = (res1.x+res1.cor(:,2)).*pixelsize;
    res1.y = (res1.y+res1.cor(:,1)).*pixelsize;
    res1.phi = z_ang;
    res1.rms = rms;
    res1.rmp = rmp;
    res1.ct = ang_ctr;
    res1.t = res1.cor(:,3);
    res1.f = ones(size(res1.x)).*ff;
    
    Res = catstruct(Res,res1,1);
    toc
end




