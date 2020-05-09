function [Pos,spL] = refine_zpos(res3,obj,zT_fit,phi0_fit,P_phi,plotflag)
P = res3.P;
Zphi = res3.phi;
coords = res3.cor;
Sx = res3.sx;
Sy = res3.sy;


pixelsize = obj.Pixelsize.*1e3;
x_fit = (P(:,1)+coords(:,1)).*pixelsize;
y_fit = (P(:,2)+coords(:,2)).*pixelsize;
z_fit = P(:,3).*1e3;
Pos = cat(2,x_fit,y_fit,z_fit);

sr_fit = Sx.^2-Sy.^2;
z_phifit = Zphi;
y = z_fit.*2*pi/zT_fit+phi0_fit;
if plotflag == 1
    figure;plot(z_fit,z_phifit,'.',z_fit,wrapToPi(y),'.')
end
zseg = round((wrapToPi(y)-y)/2/pi);
Ns = unique(zseg);

%%
pz = double(median(P_phi)*zT_fit/2/pi);
N1 = numel(z_fit);
density = zeros(N1,1);
rsr = 0.2;
tic
parfor ii = 1:N1
    mask = ((sr_fit-sr_fit(ii)).^2+(z_fit-z_fit(ii)).^2./pz./pz) < rsr^2;
    density(ii) = sum(mask);
end
toc

%%

spL = [];
for nn = 1:numel(Ns)-1
    mask1 = zseg == Ns(nn);
    pos1 = Pos(mask1,[1,2,3]);
    phi1 = z_phifit(mask1);
    sr1 = sr_fit(mask1);
    des1 = density(mask1);
    
    mask2 = zseg == Ns(nn+1);
    pos2 = Pos(mask2,[1,2,3]);
    phi2 = z_phifit(mask2);
    sr2 = sr_fit(mask2);
    des2 = density(mask2);
    
    Dxy = pdist2(pos1(:,[1,2]),pos2(:,[1,2]));
    Dz = pdist2(pos1(:,3),pos2(:,3));
    
    maskg = abs(Dz-zT_fit)<25 & abs(Dxy)<30;
    if plotflag == 1
    dipshow(maskg)
    figure;plot(Dxy(maskg),Dz(maskg),'.');axis equal
    end
    %%
    pos1_refine = pos1;
    pos2_refine = pos2;
    var=10;
    while var>0
        %%
        %close all
        [var,ind] = max(sum(maskg,1));
        ind1 = find(maskg(:,ind)==1);
        ind2 = [];
        for ii = 1:numel(ind1)
            
            indi = find(maskg(ind1(ii),:) == 1);
            ind2 = cat(2,ind2,indi);
            
        end
        ind2 = unique(ind2);
        
        for ii = 1:numel(ind2)            
            indi = find(maskg(:,ind2(ii)) == 1);
            ind1 = cat(1,ind1,indi);
        end
        ind1 = unique(ind1);
        
        avgdes1 = median(des1(ind1));
        avgdes2 = median(des2(ind2));
        avgz1 = median(pos1(ind1,3));
        avgz2 = median(pos2(ind2,3));

        avgsr1 = median(sr1(ind1));
        avgsr2 = median(sr2(ind2));
%         avgsr1 = sum(sr1(ind1).*des1(ind1)./sum(des1(ind1)));
%         avgsr2 = sum(sr2(ind2).*des2(ind2)./sum(des2(ind2)));
        sp = (avgz1-avgz2)/(avgsr1-avgsr2);
        spL = cat(1,spL,sp);
        if plotflag == 1
%             figure;plot3(pos1(ind1,1)-min(pos1(ind1,1)),pos1(ind1,2)-min(pos1(ind1,2)),pos1(ind1,3),'.',...
%                          pos2(ind2,1)-min(pos1(ind1,1)),pos2(ind2,2)-min(pos1(ind1,2)),pos2(ind2,3),'.')
            
            % figure;plot(z_fit,wrapToPi(y),'.')
            % hold on;plot(pos1(ind1,3),phi1(ind1),'.',pos2(ind2,3),phi2(ind2),'.')
            
            figure;plot(sr_fit,z_fit,'.')
            hold on;plot(sr1(ind1),pos1(ind1,3),'.',sr2(ind2),pos2(ind2,3),'.')
            text(10,double(avgz1),num2str(avgdes1,4))
            text(10,double(avgz2),num2str(avgdes2,4))
        end
        if sp>-60 || sp<-180
            if avgdes1<avgdes2
                pos1_refine(ind1,3) = pos1(ind1,3)-zT_fit;
            else
                pos2_refine(ind2,3) = pos2(ind2,3)+zT_fit;
            end
        end
        maskg1 = double(maskg);
        maskg1(ind1,ind2)=0;
        %dipshow(maskg1)
        maskg = maskg1;
    end
    if plotflag == 1
        figure;plot(sr_fit,z_fit,'.')
        hold on;plot(sr1,pos1_refine(:,3),'.',sr2,pos2_refine(:,3),'.')
        
        figure;plot(sr_fit(mask1),z_fit(mask1),'.',sr1,pos1_refine(:,3),'.')
        figure;plot(sr_fit(mask2),z_fit(mask2),'.',sr2,pos2_refine(:,3),'.')
        figure;plot(sr1,pos1_refine(:,3),'.',sr2,pos2_refine(:,3),'.')
    end
    
    Pos(mask1,3) = pos1_refine(:,3);
    Pos(mask2,3) = pos2_refine(:,3);
end

