function [estx,esty,sigmax,sigmay,zstep] = genastparam(P,stepsize,fitrange,figfolder,filename)

%% plot sigmax sigmay

xstep=(1:length(P(:,1)))'*stepsize;% nm
sigmax=P(:,5);
sigmay=P(:,6);


%% generate spline fit
stacksz=1;
sigxmean=mean(reshape(sigmax,[stacksz numel(sigmax)/stacksz]),1);

sigymean=mean(reshape(sigmay,[stacksz numel(sigmay)/stacksz]),1);

zstep=(1:length(sigxmean))*20;

sigxp=spline(zstep,sigxmean);

sigyp=spline(zstep,sigymean);


%%
start_point = [1.1 600 1200 0.1 0.1];
figure
[estx, model]=fit_ast(xstep(fitrange),sigmax(fitrange),start_point);
[~,fit_curve]=model(estx);
plot(xstep(fitrange),fit_curve,'b-','linewidth',2);
hold on
plot(xstep,sigmax,'bo');

[esty, model]=fit_ast(xstep(fitrange),sigmay(fitrange),start_point);
[~,fit_curve]=model(esty);
plot(xstep(fitrange),fit_curve,'r-','linewidth',2);
plot(xstep,sigmay,'ro');

xlabel('relative z (nm)');
ylabel('sigma');
legend('Fit. \sigma_x','Obs. \sigma_x','Fit. \sigma_y','Obs. \sigma_y','Location','North');
xlim([0, 2000]);
saveas(gcf,[figfolder,'\', filename(1:end-4), '_sigmaxy_fit_z_v1.png'],'png');

end