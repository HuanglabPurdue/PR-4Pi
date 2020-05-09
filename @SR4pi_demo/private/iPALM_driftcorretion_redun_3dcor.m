% nfs: number of frame per segment
% x,y: position in pixel
% z: positions in nm
% t: frame number
% xshift,yshift,zshift: shift in nm

function [xshift,yshift,zshift,errorf1,errorf2,errorf3] = iPALM_driftcorretion_redun_3dcor(x,y,z,pixelsize,sigma,cormask,stepsize)

z = z./pixelsize;% in pixels
G = 2;
minz = min(z);
z = z-minz+G;


errothresh = 15; % in nm
thresh = errothresh/pixelsize;
segnum = max(cormask);


if nargin > 6
    rf = 1; % redundancy factor 
else
    rf = 6;
end
plotflag = 0;
if plotflag == 1
    figure;
    ha = axes;
    hold(ha,'on');
end
shift1 = cell(1,segnum-1);
shift2 = cell(1,segnum-1);
shift3 = cell(1,segnum-1);
% parfor requires more memory
if nargin >6
    for ii = 1:segnum-1
        ii
        ind = min(ii+rf,segnum);
        maski = (cormask>=ii)&(cormask<=ind);
        xi = x(maski);
        yi = y(maski);
        zi = z(maski);
        cormi = cormask(maski);
        tic
        [dx, dy, dz] = findfineshift_3d(xi,yi,zi,cormi,sigma,stepsize);
        toc
        %     for jj = ii+1:min(ii+rf,segnum)
        %         AA(kk,ii+1:jj) = 1;
        %         kk = kk+1;
        %     end
        %     kk = kk+1;
        shift1{ii} = dx';
        shift2{ii} = dy';
        shift3{ii} = dz';
        
        
    end
    
else
    for ii = 1:segnum-1
        ii
        ind = min(ii+rf,segnum);
        maski = (cormask>=ii)&(cormask<=ind);
        xi = x(maski);
        yi = y(maski);
        zi = z(maski);
        cormi = cormask(maski);
        tic
        [dx, dy, dz] = findfineshift_3d(xi,yi,zi,cormi,sigma);
        toc
        %     for jj = ii+1:min(ii+rf,segnum)
        %         AA(kk,ii+1:jj) = 1;
        %         kk = kk+1;
        %     end
        %     kk = kk+1;
        shift1{ii} = dx';
        shift2{ii} = dy';
        shift3{ii} = dz';
        
        
    end
    
end

shift1 = cell2mat(shift1)';
shift2 = cell2mat(shift2)';
shift3 = cell2mat(shift3)';
if plotflag == 1
    drawnow;
    plot(shift1,'r.','parent',ha);
    plot(shift2,'b.','parent',ha);
    plot(shift3,'g.','parent',ha);
    legend('x','y','z')
end
clear AA
kk=2;
for ii = 1:segnum-1
    for jj = ii+1:min(ii+rf,segnum)
        AA(kk,ii+1:jj) = 1;
        kk = kk+1;
    end
    kk = kk+1;
end

R_shift1=pinv(AA)*shift1;
R_shift2=pinv(AA)*shift2;
R_shift3=pinv(AA)*shift3;

error1=sqrt(((AA*R_shift1)-shift1).^2);
error2=sqrt(((AA*R_shift2)-shift2).^2);
error3=sqrt(((AA*R_shift3)-shift3).^2);

[rshift1,errorf1]= rankshift(AA,shift1,error1,thresh);
[rshift2,errorf2]= rankshift(AA,shift2,error2,thresh);
[rshift3,errorf3]= rankshift(AA,shift3,error3,thresh);

nl = numel(rshift1);
M = tril(ones(nl,nl));

xshift = (M*rshift1).*pixelsize;
yshift = (M*rshift2).*pixelsize;
zshift = (M*rshift3).*pixelsize;


end

