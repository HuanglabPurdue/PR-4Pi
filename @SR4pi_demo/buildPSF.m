function buildPSF(obj)
% buildPSF - generate PR-4PiPSF model from phase retrieved pupil functions
%   It generates a 4PiPSF model using the PSF_4pi class given the pupil
%   functions generated from getPR() and the calibration parameters obtained
%   from caliparam(). 
%
% Output:
%   PSF4pi - 4PiPSF model
%   Phi_offset - offset of measured cavity phase to the theoretical cavity phase, it is near zero, unit: radian
%   Srcenter - shape metric for infocus astigmatism PSF. This is an initial estimation of srcenter.
%
% Output will be stored in obj
%%
namei = fields(obj.PSFcalidata);
psfm = [];
psf4dm = [];
for ii = 1:length(namei)
    psfm = cat(2,psfm,obj.PSFcalidata.(namei{ii}));
    psf4dm = cat(4,psf4dm,obj.PSFcalidata.(namei{ii}));
end
h1 = dipshow(psfm);                                     % measured 4Pi-PSF
diptruesize(h1,1000)

% generate initial 4pi PSF model
R = 128;
sigmax = obj.OTFsigma;
sigmay = obj.OTFsigma;
bxsz = size(psf4dm,1);
Num = size(psf4dm,3);
PRstruct1 = obj.PSFtop.PRobj.PRstruct;
PRstruct2 = obj.PSFbot.PRobj.PRstruct;
PRstruct1.Lambda = obj.Lambda;              % micron
PRstruct2.Lambda = obj.Lambda;              % micron
PRstruct1.SigmaX = sigmax;
PRstruct1.SigmaY = sigmay;
PRstruct2.SigmaX = sigmax;
PRstruct2.SigmaY = sigmay;
psfobj = PSF_4pi(PRstruct1);
Zpos = linspace(-1,1,Num);
psfobj.Xpos = zeros(Num,1);
psfobj.Ypos = zeros(Num,1);
psfobj.Boxsize = bxsz;
psfobj.Pixelsize = obj.Pixelsize;           % micron
psfobj.PSFsize = R;
psfobj.nMed = obj.RI_med;
psfobj.Phasediff = obj.PhaseparamM0.phip-obj.PhaseparamM0.phis;
psfobj.Iratio = 0.7;
psfobj.Phi0 = 0;
psfobj.Zpos = Zpos;
psfobj.Zoffset = 0;
psfobj.ModulationDepth = 0.8;

psfobj.gen2Pupil(PRstruct1,PRstruct2);
psfobj.genPupil_4pi('noIMMaber');
psfobj.genPSF_4pi_md();

label = {'p1','s2','p2','s1'};
psfr = [];
psf4dr = [];
for nn = 1:4
    psfobj.PSFs = psfobj.PSF4pi.(label{nn});
    psfobj.scalePSF();
    psfI = psfobj.ScaledPSFs;
    psfr = cat(2,psfr,psfI);
    psf4dr = cat(4,psf4dr,psfI);
end
h = dipshow(psfr);                                      % 4Pi-PSF model
diptruesize(h,1000)

obj.PSF4pi = struct('psfobj',psfobj,'PRstruct1',PRstruct1,'PRstruct2',PRstruct2);

%% find phi offset
Num = 1;
zpos = 0;
xpos = zeros(Num,1);
ypos = zeros(Num,1);
I = 1000;
bg = 1;
P1 = cat(2,xpos,ypos,zpos,I.*ones(Num,4),bg.*ones(Num,4));
N0 = 100;
phi0 = linspace(0,2*pi,N0);
w = [1,1,1,1,1];
psf4d = [];
for ii = 1:N0
    psfobj.Phi0 = phi0(ii);
    [psf] = genpsf4pi_real(P1,w,psfobj,PRstruct1,PRstruct2);
    psf4d = cat(3,psf4d,psf);
end

P2 = GPUgaussMLEv2(single(sum(psf4d,4)),1.4,200,4);
[angc] = getphase(P2,psf4d,0);
sx = P2(:,5);
sy = P2(:,6);
src = sx.^2 - sy.^2;

obj.Phi_offset = mean(phi0-unwrap(angc)');
obj.Srcenter = src(1);                      % initial srcenter, this might not be accurate


end