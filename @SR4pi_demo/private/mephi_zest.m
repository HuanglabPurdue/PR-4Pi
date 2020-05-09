function [zest_f,zerr_f, mephimask]=mephi_zest(mephipp,uwmephi,mtc_seg,zang_seg,centermtc,freq)
% center phi_0 can not be obtained by option 'cofocal' which will perform a
% estimate of the phi_0 based on the zang_seg data which is more
% accurate than mephipp.

% center phi_0 is the distribution center of the mephi plot
% phi_0=ppval(mephipp,(max(uwmephi(:,1))+min(uwmephi(:,1)))/2); %phi_mid can be set
if isempty(centermtc)
    centermtc=0;
    phi_0=ppval(mephipp,centermtc);
    %phi_0=wrapToPi(phi_0);
elseif strcmp(class(centermtc),'double')==1
    tmpmask=mtc_seg>(centermtc-0.02)&mtc_seg<(centermtc+0.02);
    tmp=zang_seg(tmpmask);
    [tmphis xout]=hist(tmp,20);
    inis=[max(tmphis(:)) mean(tmp(:)) 0.5];
    [c]=fit1dgaussian(tmphis,xout,inis);
    phi_0=c(2);
elseif strcmp(centermtc,'mid')
    [nout xout]=hist(mtc_seg,20);
    [I ind]=max(nout);
    centermtc2=xout(ind);
    phi_0=ppval(mephipp,centermtc2);   
elseif strcmp(centermtc,'cofocal')
    tmpmask=mtc_seg>-0.02&mtc_seg<0.02;
    tmp=zang_seg(tmpmask);
    [tmphis xout]=hist(tmp,20);
    inis=[max(tmphis(:)) mean(tmp(:)) 0.5];
    [c]=fit1dgaussian(tmphis,xout,inis);
    phi_0=c(2);
end

if phi_0>=6*pi||phi_0<=-6*pi
    display(['The phi_0 estimation did not converge, setting this value to zero;']);
    phi_0=0;
end


% phi_0=ppval(mephipp,centermtc); %phi_mid set to 0 for segment align

display(['phi_0 is ' num2str(phi_0)]);

tmp=-10*2*pi:2*pi:10*2*pi;
pispace=repmat(tmp,[size(zang_seg,1) 1]);
angle1space=repmat(zang_seg,[1 size(pispace,2)]);
angle1span=(angle1space+pispace);%data
phispan=ppval(mephipp,mtc_seg);
phispan=repmat(phispan(:),[1 size(pispace,2)]);%model
phidiff1=abs(phispan-angle1span);

phidiff=phidiff1;
[C, I]=min(phidiff,[],2);
mindiffspan=repmat(C,[1 size(pispace,2)]);
mask=(mindiffspan==phidiff);
uwangle=sum((double(mask).*angle1span),2);
angle_error=C;

norm_uwangle=(uwangle-phi_0);

zest=norm_uwangle./freq;
zerr=angle_error./freq;
%filtering
mephimask=mtc_seg<max(uwmephi(:,1))&mtc_seg>min(uwmephi(:,1));
zerrmask = zerr<80; % nm
mask = mephimask&zerrmask;
zest_f = zest;
zerr_f = zerr;

h = figure;
h.Position = [1162, 598, 560, 420];
plot(uwmephi(:,1),uwmephi(:,2),'o')
hold on
plot(mtc_seg,phispan(:,1),'g.','markersize',2)
plot(mtc_seg,norm_uwangle,'.','markersize',1)
plot(mtc_seg(mask),norm_uwangle(mask),'c.','markersize',1)
legend('mephi','interpolation fit','data')
xlabel('sigma metric')
ylabel('unwrapped phase')

