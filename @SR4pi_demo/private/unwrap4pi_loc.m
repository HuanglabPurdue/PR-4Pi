function [psf4dt_peak,x0_peak,psft_peak,cor_peak,x_phifit,y_phifit,z_phifit,maskp] = unwrap4pi_loc(xp_sr,indm,p_phi,phi_thresh,astcal,sr_data,phi_data,psf4dt_data,x0_data,cor_data,plotflag)
Np = numel(xp_sr);
phi_peak = zeros(numel(sr_data),Np);
for ii = 1:Np
    phi_peak(:,ii) = polyval([p_phi,0],sr_data)-p_phi*xp_sr(ii);
    
end
if plotflag == 1
    figure;plot(sr_data,phi_data,'.')
    hold on
    plot(sr_data,phi_peak)
end

phidataL = repmat(phi_data,1,Np);
phidiff = phidataL - phi_peak;
mask = abs(phidiff)<phi_thresh;
maskp = sum(mask,2)==1;
phidata_peak = phi_data(maskp);
srdata_peak = sr_data(maskp);

if plotflag == 1
    figure;plot(sr_data,phi_data,'.')
    hold on
    plot(srdata_peak,phidata_peak,'.')
end

indT = [1:Np]-indm;
indL = repmat(indT,numel(sr_data),1);
ind_data = sum(indL.*double(mask),2);
ind_peak = ind_data(maskp);
phidata_unwrap = phidata_peak-2*pi.*ind_peak;
if plotflag == 1
    figure;plot(sr_data,phi_data,'.')
    hold on
    plot(srdata_peak,phidata_unwrap,'.')
end
psf4dt_peak = psf4dt_data(:,:,maskp,:);
x0_peak = x0_data(maskp,:);
psft_peak = cat(2,psf4dt_peak(:,:,:,1),psf4dt_peak(:,:,:,2),psf4dt_peak(:,:,:,3),psf4dt_peak(:,:,:,4));
cor_peak = cor_data(maskp,:);

xdata_peak = x0_peak(:,1)+cor_peak(:,1);
ydata_peak = x0_peak(:,2)+cor_peak(:,2);
x_phifit = xdata_peak.*129;
y_phifit = ydata_peak.*129;
z_phifit = phidata_unwrap.*astcal.zT./2./pi;
if plotflag == 1
    figure;
    scatter3(x_phifit,y_phifit,z_phifit,20,z_phifit,'.')
    axis equal
end


