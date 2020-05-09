function [p1,astparam] = genzestparam(estx,esty,sigmax,sigmay,zstep,fitrange,figfolder,filename)

Ax = estx(4);
Ay = esty(4);
Bx = estx(5);
By = esty(5);
PSFsigmax = estx(1);
PSFsigmay = esty(1);
gamma = (estx(2) - esty(2))/2/1e3;
d = (estx(3) + esty(3))/2/1e3;
Sx = sigmax;
Sy = sigmay;
z = zstep./1e3;
z0 = 1;%um
startpoint=double([Ax,Bx,Ay,By,gamma,d,PSFsigmax,PSFsigmay,z0]);

astparam=fminsearch(@(x) AstCalibr(x,Sx(fitrange),Sy(fitrange),z(fitrange),1),startpoint);
[sse,Sx1,Sy1]=AstCalibr(astparam,Sx,Sy,z,1);
figure('position',[200,300,500,400],'color',[1,1,1])
plot(z,Sx,'r.',z,Sy,'b.')
hold on
plot(z,Sx1,'r-',z,Sy1,'b-')
legend('found \sigma_x','found \sigma_y','calibration curve (\sigma_x)','calibration curve (\sigma_y)')
xlabel('z (\mum)','fontsize',12)
ylabel('\sigma_x, \sigma_y, (pixel)','fontsize',12)
set(gca,'fontsize',12)
saveas(gcf,[figfolder,'\', filename(1:end-6), '_sigmaxy_fit_z_v2.png'],'png');


Sr=Sx.^2-Sy.^2;
Srfit = Sr(fitrange);
zfit = z(fitrange);
p1=polyfit(Srfit,zfit,3);
f1=polyval(p1,Srfit);
figure('position',[200,200,700,600]);
plot(Sr,z,'o','markersize',10);
hold on
plot(Srfit,f1,'r.','linewidth',2,'markersize',10);
title('Calibration curve for initial z estimation','fontsize',12)
xlabel('Sigma R (pixels)','fontsize',12)
ylabel('z (\mum)','fontsize',12)
set(gca,'fontsize',12)
saveas(gcf,[figfolder,'\', filename(1:end-6), '_z_estimation.png'],'png');
end