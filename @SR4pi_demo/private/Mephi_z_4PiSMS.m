function [zest, zerr, zmask]=Mephi_z_4PiSMS(res,fnum,stopval,centermtc,freq)
zangresult = res.phi;
mtcresult = res.mtc;
maxf = max(res.t)+1;
tresult = res.t+(res.f-1)*maxf;


% setup
msz = 256;
srange = [0.2 0.6];
maxangle = 0.2;
sigma = 7; % from 10, for Spermatocytes to 5 for cy3b

segnum = ceil((max(tresult))/fnum);
zest = [];
zerr = [];
zmask = [];
for ii = 1:segnum
    tic
    st=(ii-1)*fnum;
    if ii==segnum
        ed=max(tresult);
    else
        ed=(ii)*fnum-1;
    end
    
    maskt=tresult>=st&tresult<=ed;
    mtc_seg=mtcresult(maskt);
    zang_seg=zangresult(maskt);
    [dmap]=build_dmap(mtc_seg,zang_seg,msz,sigma);
    h = dipshow(dmap);
    h.Position = [459   797   256   256];
    [mephi_ini,cpeak]=find_mephi_ini(dmap,srange,msz);
    stopval=max(stopval,cpeak./25);
    [mephi,uwmephi,mephipp]=find_mephi(dmap,mephi_ini,maxangle,srange,stopval);
    [zest_f,zerr_f,mephimask]=mephi_zest(mephipp,uwmephi,mtc_seg,zang_seg,centermtc,freq);

    zest=cat(1,zest,zest_f(:));
    zerr=cat(1,zerr,zerr_f(:));
    zmask=cat(1,zmask,mephimask(:));
    toc
    
    h = figure;
    h.Position = [1131, 127, 400, 350];
    scatter(mtc_seg,zang_seg,2)
    hold on
    scatter(mephi(:,1),mephi(:,2),25,'linewidth',1.5)
    xlim([-pi pi]);
    ylim([-1*pi 1*pi]);
    xlabel('sigma metric')
    ylabel('unwrapped phase')

    pause(1)
    close all
end

  