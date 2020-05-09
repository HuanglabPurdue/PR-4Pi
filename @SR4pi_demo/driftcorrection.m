function driftcorrection(obj,type)
% driftcorrection - post drift correction of localization results
%
% Input:
%   type - 'PR': drift correction for localization results from PR-4Pi algorithm
%          'CT': drift correction for localization results from contrast algorithm
%
% Output:
%   ResultPR.Final - drift correction results for PR-4Pi localizations
%       res - result structure, see mlefit() for details
%       shift - x,y,z shifts of each data segment relative to the first data segment within one optical section
%       shift_stack - x,y,z shifts of each optical section relative to the first optical section, only output for multiple optical section imaging
%       shift_stack_ast - similar as above but estimated based on localizations with astigmatism method, which tend to give a smaller shift comparing with 4Pi localization results
%   ResultCT.Final - drift correction results for PR-4Pi localizations
%       res - result structure, see initfit() and unwrap() for details
%       shift - x,y,z shifts of each data segment relative to the first data segment within one optical section
%       shift_stack - x,y,z shifts of each optical section relative to the first optical section, only output for multiple optical section imaging
%       shift_stack_ast - similar as above but estimated based on localizations with astigmatism method, which tend to give a smaller shift comparing with 4Pi localization results
%
% Output will be stored in obj and saved in file with key word '_dc_single'
% or '_dc_stack' under obj.Resultpath\PR fit\ or obj.Resultpath\CT fit\
%%
switch type
    case 'PR'
        figfolder=[obj.Resultpath,'\PR fit\'];
    case 'CT'
        figfolder=[obj.Resultpath,'\CT fit\'];
end

if ~exist(figfolder,'dir')
mkdir(figfolder)
end
%%
pixelsize = obj.Pixelsize.*1e3;
n_med = obj.RI_med;
n_imm = obj.RI_imm;

switch type
    case 'PR'
        res = obj.ResultPR.MLE.res;
        res.cycN = res.cycN-min(res.cycN)+1;
        res.stepN = res.stepN-min(res.stepN)+1;
        maskzsr = obj.ResultPR.MLE.maskrm;
        res_fit0 = rmfield(res,{'x0','z_guess','crlb','transx','transy'});
        res_fit = applymask(res_fit0,maskzsr);
        llthresh = 2e3;
    case 'CT'
        res_fit = obj.ResultCT.unwrap;
        res_fit.cycN = res_fit.cycN-min(res_fit.cycN)+1;
        res_fit.stepN = res_fit.stepN-min(res_fit.stepN)+1;
        llthresh = 600;
end

Ns = max(res_fit.stepN);
xLim = [30,140].*pixelsize;
yLim = [30,155].*pixelsize;
zLim = [-600, 600];
dcflag = 1;
if dcflag == 1
    shift_fit = cell(1,Ns);
end
for ii = 1:Ns
    mask_step = res_fit.stepN==ii;
    res_step = applymask(res_fit,mask_step);
    
    maxf = max(res_step.t)+1;
    frames = res_step.t+(res_step.cycN-1).*maxf;
    masksubx = res_step.x>xLim(1) & res_step.x<xLim(2);
    masksuby = res_step.y>yLim(1) & res_step.y<yLim(2);
    maskz = res_step.z>zLim(1) & res_step.z<zLim(2);
    maskLL = res_step.LL<llthresh;
    masksub = masksubx & masksuby & maskLL & maskz;
    fnum = 1*maxf;
    
    resi = applymask(res_step,masksub);
    
    sigma = 2;
    [cormask] = gencormask(frames,fnum);
    [cormasks] = gencormask(frames(masksub),fnum);
    
    if dcflag == 1
        [dx1,dy1,dz1] = iPALM_driftcorretion_redun_3dcor((resi.x-xLim(1))./pixelsize,(resi.y-yLim(1))./pixelsize,resi.z,pixelsize,sigma,cormasks);
        figure;plot([dx1,dy1,dz1],'o-')
        legend('x','y','z')
        xlabel('batch number')
        ylabel('drift (nm)')
        shift_fit{ii} = [dx1,dy1,dz1];
    end
    [xout,yout,zout] = shiftcoords(res_step.x,res_step.y,res_step.z,cormask,-1.*shift_fit{ii});
    [~,~,zout_ast] = shiftcoords([],[],res_step.zast,cormask,-1.*shift_fit{ii});
    
    res_fit.x(mask_step) = xout;
    res_fit.y(mask_step) = yout;
    res_fit.z(mask_step) = zout;
    res_fit.zast(mask_step) = zout_ast;
end
switch type
    case 'PR'
        obj.ResultPR.Final.res = res_fit;
        obj.ResultPR.Final.shift = shift_fit;
    case 'CT'
        obj.ResultCT.Final.res = res_fit;
        obj.ResultCT.Final.shift = shift_fit;
end
res = res_fit;
shift = shift_fit;
save([figfolder,obj.Savename,'_dc_single.mat'],'res','shift')
if Ns>1
    %% stack alignment
    zLim = [-800, 800];
    masksubx = res_fit.x>xLim(1) & res_fit.x<xLim(2);
    masksuby = res_fit.y>yLim(1) & res_fit.y<yLim(2);
    maskz = res_fit.z>zLim(1) & res_fit.z<zLim(2);
    maskzast = res_fit.zast>zLim(1) & res_fit.zast<zLim(2);
    maskLL = res_fit.LL<llthresh;
    masksub = masksubx & masksuby & maskLL & maskz & maskzast;
    
    resi = applymask(res_fit,masksub);
    cormask = res_fit.stepN;   % cormask should be 1 based
    cormasks = res_fit.stepN(masksub);
    stepsize = obj.Stepsize.*1e3*n_med/n_imm/pixelsize; % in pixel, step size should be changed for each data
    
    if dcflag == 1
        [dx1,dy1,dz1] = iPALM_driftcorretion_redun_3dcor((resi.x-xLim(1))./pixelsize,(resi.y-yLim(1))./pixelsize,resi.z,...
            pixelsize,sigma,cormasks,stepsize);
        shift_stack = [dx1,dy1,dz1];

        [dx2,dy2,dz2] = iPALM_driftcorretion_redun_3dcor((resi.x-xLim(1))./pixelsize,(resi.y-yLim(1))./pixelsize,resi.zast,...
            pixelsize,sigma,cormasks,stepsize);
        shift_stack_ast = [dx2,dy2,dz2];
        figure;plot(shift_stack,'o-');hold on;
        plot(shift_stack_ast,'o-')
        legend('x','y','z','x ast.','y ast.','z ast.')
        xlabel('optical section number')
        ylabel('drift (nm)')

    end
    [xout,yout,zout] = shiftcoords(res_fit.x,res_fit.y,res_fit.z,cormask,-1.*shift_stack);
    [~,~,zout_ast] = shiftcoords([],[],res_fit.zast,cormask,-1.*shift_stack_ast);
    
    res_fit.x = xout;
    res_fit.y = yout;
    res_fit.z = zout;
    res_fit.zast = zout_ast;
    switch type
        case 'PR'
            obj.ResultPR.Final.res = res_fit;
            obj.ResultPR.Final.shift_stack = shift_stack;
            obj.ResultPR.Final.shift_stack_ast = shift_stack_ast;
        case 'CT'
            obj.ResultCT.Final.res = res_fit;
            obj.ResultCT.Final.shift_stack = shift_stack;
            obj.ResultCT.Final.shift_stack_ast = shift_stack_ast;
    end
    res = res_fit;
    save([figfolder,obj.Savename,'_dc_stack.mat'],'res','shift_stack','shift_stack_ast')
end

end