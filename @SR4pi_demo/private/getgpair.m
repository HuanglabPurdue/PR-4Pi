function [inifit,gpair,z_guess] = getgpair(res,gzone,pixelsize,plotflag)

PcudaM = res.P;
errM = res.err;
cor_peak = res.cor;
z_ast_peak = res.zast;
inifit.x = (PcudaM(:,1)+cor_peak(:,1)).*pixelsize.*1e3;
inifit.y = (PcudaM(:,2)+cor_peak(:,2)).*pixelsize.*1e3;
inifit.z = PcudaM(:,3).*1e3;
inifit.zast = z_ast_peak.*1e3;
if size(PcudaM,2)>5
    inifit.I = sum(PcudaM(:,4:7),2);
else
    inifit.I = PcudaM(:,4).*4;
end
    
inifit.LL = errM(:,2);
Nfit = numel(inifit.x);
z_guess = inifit.z;

if plotflag == 1
    disx = pdist(cat(2,inifit.x,inifit.y),'euclidean');
    disz = pdist(cat(2,inifit.z,zeros(Nfit,1)),'euclidean');
    h = figure;
    ha = axes;
    hold on
    plot(disx,disz,'b.','markersize',1);axis equal;
    ha.XLim = [0,100];
    ha.YLim = [200,400];
end

if Nfit<10000
    [gpair,~] = findgpair(inifit,gzone,plotflag);
    [z_guess] = refineini4pi_z(inifit,gzone,gpair);
else
    gpair = [];
    Ns = 2;
    xlim0 = linspace(min(inifit.x), max(inifit.x),Ns+1);
    ylim0 = linspace(min(inifit.y), max(inifit.y),Ns+1);
    for ii = 1:Ns
        for jj = 1:Ns
            masksubx1 = inifit.x>=xlim0(jj) & inifit.x<xlim0(jj+1);
            masksuby1 = inifit.y>=ylim0(ii) & inifit.y<ylim0(ii+1);
            masksub = masksubx1 & masksuby1;
            res_fitsub = applymask(inifit,masksub);
            if numel(res_fitsub.x) > 2
                [gpairi,~] = findgpair(res_fitsub,gzone,plotflag);
                [z_guess1] = refineini4pi_z(res_fitsub,gzone,gpairi);
                gpair = cat(1,gpair,gpairi);
                z_guess(masksub) = z_guess1;
            end
        end
    end
end
