        
function [P1,err1,cg1,crlb1,std1] = loc4pi_cuda_sCMOS_constI(data,gainR,obj,res,rt,plotflag)
sobj = obj.PSFlib;
lambda = obj.Dampingfactor;
iterateN = obj.IterationNum;
boxsize = obj.Boxsize;
cor = res.cor;
x0 = res.x0;
trans_dx = res.transx;
trans_dy = res.transy;
z_guess = res.z_guess;
N = size(data,3);
maxN = 10000;
v0 = [1:maxN:N];
Ns = numel(v0);
if v0(Ns)<N
    v0(Ns+1) = N;
end
v0(end) = v0(end)+1;
xtmp = x0(:,1)+boxsize/2;
ytmp = x0(:,2)+boxsize/2;
if size(x0,2)>5
    x0 = cat(2,xtmp,ytmp,z_guess,x0(:,4:11));
else
    x0 = cat(2,xtmp,ytmp,z_guess,x0(:,4:5));
end
tic
P1 = [];
err1 = [];
cg1 = [];
crlb1 = [];
std1 = [];
tic
for ii = 1:numel(v0)-1
    vi = v0(ii):v0(ii+1)-1;
    [PcudaM,cgM,crlbM,errM,stdM] = genloc(data(:,:,vi,:),gainR(:,:,vi,:),cor(vi,:),x0(vi,:),sobj,iterateN,lambda,trans_dx(vi,:),trans_dy(vi,:),rt,plotflag);
    P1 = cat(1,P1,PcudaM);
    err1 = cat(1,err1,errM);
    cg1 = cat(1,cg1,cgM);
    crlb1 = cat(1,crlb1,crlbM);
    std1 = cat(1,std1,stdM);
end
time1 = toc;
disp(['Elapsed time for MLE (cubic spline): ',num2str(time1),' s'])
end


function [PcudaM,cgM,crlbM,errM,stdM] = genloc(data,gainR,cor,x0,sobj,iterateN,lambda,trans_dx,trans_dy,rt,plotflag)
data_cuda = single(data(:));
gainR_cuda = single(gainR(:));
bxsz = size(data,1);
coords_cuda = single(cor(:));
Nfit = numel(x0(:,1));
if size(x0,2)>5
    x01 = cat(2,x0(:,1:3),mean(x0(:,4:7),2),mean(x0(:,8:11),2));
else
    x01=x0;
end
x0i = single(reshape(x01',Nfit*5,1));


[P_cuda,CG,crlb,err,PSF_cuda] = cuda4pi_loc_spline_consI(data_cuda,coords_cuda,single(sobj.st.F),...
    sobj.dx,sobj.dz,sobj.startx,sobj.starty,sobj.startz,iterateN,Nfit,lambda,x0i,gainR_cuda,...
    single(sobj.st.Fx),single(sobj.st.Fy),single(sobj.st.Fz),...
    single(sobj.st.Fxy),single(sobj.st.Fxz),single(sobj.st.Fyz),single(sobj.st.Fxyz),...
    single(trans_dx(:)),single(trans_dy(:)));

clear cuda4pi_loc_spline_consI

crlbM = reshape(crlb,5,Nfit)';
stdM = sqrt(crlbM);
errM = reshape(err,2,Nfit)';
cgM = reshape(CG,5,Nfit)';
PcudaM = reshape(P_cuda,5,Nfit)';
tmp_psf = permute(reshape(PSF_cuda,[bxsz,bxsz,Nfit,4]),[1,2,3,4]);
tmp_data = permute(reshape(data_cuda,[bxsz,bxsz,Nfit,4]),[1,2,3,4]);

if plotflag == 1
    psf_cuda = [];
    Data_cuda = [];
    for ss = 1:4
        psf_cuda = cat(2,psf_cuda,tmp_psf(:,:,:,ss));
        Data_cuda = cat(2,Data_cuda,tmp_data(:,:,:,ss));
    end
    ov = joinchannels('RGB',Data_cuda,psf_cuda);
    dipshow(ov)
    dipshow(cat(1,Data_cuda,psf_cuda))
end


% im_fit = zeros(168,168,2000,4);
% im_data = zeros(168,168,2000,4);
% s1 = floor(bxsz/2);
% s2 = ceil(bxsz/2);
% for ii = 1:Nfit
%     for ss = 1:4
%         im_fit(cor(ii,2)-s1+1:cor(ii,2)+s2,cor(ii,1)-s1+1:cor(ii,1)+s2,cor(ii,3)+1,ss) = tmp_psf(:,:,ii,ss)+im_fit(cor(ii,2)-s1+1:cor(ii,2)+s2,cor(ii,1)-s1+1:cor(ii,1)+s2,cor(ii,3)+1,ss);
%         im_data(cor(ii,2)-s1+1:cor(ii,2)+s2,cor(ii,1)-s1+1:cor(ii,1)+s2,cor(ii,3)+1,ss) = tmp_data(:,:,ii,ss)+im_data(cor(ii,2)-s1+1:cor(ii,2)+s2,cor(ii,1)-s1+1:cor(ii,1)+s2,cor(ii,3)+1,ss);
%     end
% end

if plotflag == 1
    PSFsigma = 1.15;
    iterations = 100;
    fittype = 4;
    datasum = squeeze(sum(tmp_psf,4));
    datasum(datasum<=0) = 1e-6;
    [P, CRLB, LL] = GPUgaussMLEv2(single(datasum),PSFsigma,iterations,fittype);
    
    sx1 = P(:,5);
    sy1 = P(:,6);%
    sr1 = sx1.^2 - sy1.^2;
    
    datasum = squeeze(sum(tmp_data,4));
    datasum(datasum<=0) = 1e-6;
    [P, CRLB, LL] = GPUgaussMLEv2(single(datasum),PSFsigma,iterations,fittype);
    
    sx2 = P(:,5);
    sy2 = P(:,6);
    sr2 = sx2.^2 - sy2.^2;
    
    [val,ind] = min(abs(PcudaM(:,3)));
    
    h = figure;
    h.Position = [100,200,500,400];
    ha = axes;hold on;
    plot(PcudaM(:,3),sr1,'.',PcudaM(:,3),sr2,'.')
    plot(PcudaM(ind,3),sr2(ind),'mo')
    text(0.1,0.2,['sr center: ',num2str(sr2(ind),1)])
    ha.XLabel.String = 'z postion';
    ha.YLabel.String = 'sr';
    ha.YLim = [-15,8];
    axis tight
    legend('fit','data')
    
end


end


