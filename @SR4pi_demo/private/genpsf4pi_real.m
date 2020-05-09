function [psf4pi,psf4pi0] = genpsf4pi_real(x,w,obj,PRstruct1,PRstruct2,stagepos,chamberH,abertype)
R = obj.Boxsize;
obj.Xpos = x(:,1).*w(1);
obj.Ypos = x(:,2).*w(2);
N = numel(x(:,1));
% I = x(:,4);
% tmp = zeros(1,1,N);
% tmp(1,1,:) = I;
% IL = repmat(tmp,[R,R,1]);
obj.precomputeParam();
obj.gen2Pupil(PRstruct1,PRstruct2);
if nargin>5
    obj.StagePosUp = chamberH-stagepos;
    obj.StagePosDown = stagepos;
    %obj.StagePosDown = 1;
    obj.Zpos = x(:,3).*w(3);
    obj.genPupil_4pi(abertype)
else
    obj.Zpos = x(:,3).*w(3);
    obj.genPupil_4pi('noIMMaber')
end
obj.genPSF_4pi_md();
label = {'p1','s2','p2','s1'};
psf4pi = zeros(R,R,N,4);
psf4pi0 = zeros(R,R,N,4);
for nn = 1:4
    if strcmp(label{nn},'s1')
        obj.PlaneDis = 0.0; % um
    else
        obj.PlaneDis = 0;
    end
    %obj.genPSF(obj.Pupil4pi.(label{nn}));
    obj.PSFs = obj.PSF4pi.(label{nn});
    obj.scalePSF();
    psfI = obj.ScaledPSFs;
    %psf = psfI./4;
    psf = psfI;
    
    I = x(:,nn+3);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = I;
    IL = repmat(tmp,[R,R,1]).*w(4);

    bg = x(:,nn+7);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = bg;
    bgL = repmat(tmp,[R,R,1]).*w(5);
    
    psf4pi(:,:,:,nn) = psf.*IL+bgL;
    psf4pi0(:,:,:,nn) = obj.PSFs./4.*IL+bgL;
end

end