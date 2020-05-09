function [z_guess,mcc_max] = genini4pi_z_cuda(data,sobj)
N = size(data,3);
maxN = 1000;
v0 = [1:maxN:N];
Ns = numel(v0);
if v0(Ns)<N
    v0(Ns+1) = N;
end
v0(end) = v0(end)+1;
tic
z_guess = [];
mcc_max = [];
for ii = 1:numel(v0)-1
    [z0,mcc] = genini_z(data(:,:,v0(ii):v0(ii+1)-1,:),sobj);
    z_guess = cat(1,z_guess,z0);
    mcc_max = cat(1,mcc_max,mcc);
end
time1 = toc;
disp(['Elapsed time for z initial guess: ',num2str(time1),' s'])
end

function [z_guess,mcc_val] = genini_z(data,sobj)
Num = 61;
boxsize = size(data,1);
Nfit = size(data,3);
Npsf = Num*Nfit;
zs = linspace(-1,1,Num);
z = repmat(zs,Nfit,1)';
z = z(:);
x = zeros(Npsf,1)+boxsize/2;
y = zeros(Npsf,1)+boxsize/2;
I = 1;
bg = 0;
coords_cuda = cat(2,ones(Npsf,2),zeros(Npsf,1))';
coords_cuda = single(coords_cuda(:));

x0 = cat(2,x,y,z,ones(Npsf,4).*I,ones(Npsf,4).*bg);
x0i = single(reshape(x0',Npsf*11,1));
data_cuda = zeros(boxsize,boxsize,Npsf*4);
data_cuda = single(data_cuda(:));
varsub = single(zeros(size(data_cuda)));
gainsub = single(ones(size(data_cuda)));
gainRsub = varsub./gainsub./gainsub;

% simulate data
iterateN = 0;
[~,~,~,~,PSF_cuda] = cuda4pi_loc_linear(data_cuda,coords_cuda,sobj.samplepsf_cuda,...
                                sobj.dx,sobj.dz,sobj.startx,sobj.starty,sobj.startz,iterateN,Npsf,0,x0i,gainRsub);

clear cuda4pi_loc_stream

psf_model = reshape(PSF_cuda,[boxsize,boxsize,Npsf,4]);

z_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
parfor nn = 1:Nfit
    mcc_max = 0;
    for ii = 1:Num
        mcc = 0;
        for ss = 1:4
            ref = psf_model(:,:,(nn-1)*Num+ii,ss);
            img = data(:,:,nn,ss);
            mcc = mcc + cc2(ref,img)/4;
        end
        if mcc>mcc_max
            ind = ii;
            mcc_max = mcc;
        end
    end
    z_guess(nn) = zs(ind);
    mcc_val(nn) = mcc_max;
end
end
%figure;plot(z_guess)

