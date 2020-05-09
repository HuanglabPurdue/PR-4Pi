function [z_guess] = refineini4pi_z(inifit,gzone,gpair)
z_guess = inifit.z;
count = 0;
for id = 1:length(gpair)
    
    z1 = inifit.z(gpair{id});
    zmin = min(z1);
    zmax = max(z1);
    Ng = round((zmax-zmin)/gzone.zlim(2))+1;
    ind_tmp = cell(Ng,1);
    ill = zeros(Ng,1);
    propz = zeros(Ng,1);
    sdz = 40;
    hs = histogram(z1,Ng,'BinLimits',[zmin-10,zmax+10]);
    for kk = 1:Ng
        ind_tmp{kk} = find(z1>=hs.BinEdges(kk) & z1<hs.BinEdges(kk+1));
        cczast = mean(inifit.zast(gpair{id}(ind_tmp{kk})));
        ccz = mean(inifit.z(gpair{id}(ind_tmp{kk})));
        ill(kk) = max(inifit.I(gpair{id}(ind_tmp{kk}))./inifit.LL(gpair{id}(ind_tmp{kk})));
        propz(kk) = exp(-(cczast-ccz)^2/2/sdz/sdz);
    end
    [~,indm] = max(ill);
    gpairN = gpair{id}(ind_tmp{indm});
    sdx = min([max([10,std(inifit.x(gpairN))]),30]);
    sdy = min([max([10,std(inifit.y(gpairN))]),30]);
    sdz = min([max([40,std(inifit.zast(gpairN))]),100]);
%     if sdx == 0
%         sdx = 10;
%         sdy = 10;
%         sdz = 40;
%     end
    mux = median(inifit.x(gpairN));
    muy = median(inifit.y(gpairN));
    muzast = median(inifit.zast(gpairN));
    muz = median(inifit.z(gpairN));
    
    ind_tmp1 = ind_tmp;
    ind_tmp1(indm) = [];
    
    for gg = 1:length(ind_tmp1)
        gpairT = gpair{id}(ind_tmp1{gg});
        xt = mean(inifit.x(gpairT));
        yt = mean(inifit.y(gpairT));
        zast = mean(inifit.zast(gpairT));
        zt = mean(inifit.z(gpairT));
        prop = exp(-(xt-mux)^2/2/sdx/sdx)*exp(-(yt-muy)^2/2/sdy/sdy)*exp(-(zast-muzast)^2/2/sdz/sdz);
        propz1 = exp(-(zast-muzast)^2/2/sdz/sdz);
        propz2 = exp(-(zast-zt)^2/2/sdz/sdz);
        if (prop>1e-3)||(propz1>0.1)
            % it belongs to the same cluster
            z_guess(gpairT) = muz;
            count = count+1;
            %figure;plot(inifit.x(gpairT),inifit.zast(gpairT),'o',inifit.x(gpairT),z_guess(gpairT),'o');
        elseif (sign(zt-muzast)~=sign(zast-muzast))&&(propz2<1e-3)    
            z_guess(gpairT) = muz + sign(zast-muz)*gzone.zT;
        end
        
        
    end
end