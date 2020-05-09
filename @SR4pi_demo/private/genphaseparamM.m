function [phi_s,phi_p,ang_coeff,contrast,angfit,ang] = genphaseparamM(P,subims,M,plotflag)
%% 
xstep=(1:length(P(:,1)))'*20;
xco=P(:,2);
yco=P(:,1);

%% get rms rmp
xf=xco;
yf=yco;
[rms, rmp]=iPALMast_findmom_givenM(subims,xf,yf,M);

%% calculate reduced moment
rm1=rms./sqrt(rms.^2+rmp.^2);
rm2=rmp./sqrt(rms.^2+rmp.^2);
if plotflag == 1
    figure;
    plot(rm1,'r');
    hold on;
    plot(rm2,'b');
    % plot(-rm1,'--r');
    legend('rm1','rm2');
    xlabel('step number')
    ylabel('normalized intensity contrast')
end
%% estimate the phase delay
fitrange = [20:90];
A1 = 0.2;
w = 0.5;
phi1 = 1.2;
b1 = 0.1;
A2 = A1;
phi2 = phi1+1.4;
b2 = b1;
%
rmsin=rms(fitrange);
rmpin=rmp(fitrange);
x=1:numel(rmsin);
start_point=[A1 w phi1 b1 A2 phi2 b2];
[est, model]=iPALM_find_delay(x,rmsin,rmpin,start_point);
contrast = (abs(est(1))+abs(est(5)))/2;
[~,fitrms,fitrmp]=feval(model,est);

if plotflag == 1
    figure
    plot(x,rmsin,'r');
    hold on
    plot(x,rmpin,'b');
    plot(x,fitrms,'--r');
    hold on
    plot(x,fitrmp,'--b');
    legend('rms','rmp','fit rms','fit rmsp')
    xlabel('step number')
    ylabel('intensity contrast')
end

phi_s = est(3);
phi_p = est(6);
dphi = phi_s-phi_p;


%% get angle from rms and rmp

c=rm1-1i*rm2;
ang=angle(c);   % change to MLE
ang2=[];
for ii=1:numel(rms)
[ang2(ii,:)]=iPALM_est_angle(rms(ii),rmp(ii),phi_s,phi_p,ang(ii)-phi_s);
end

if plotflag == 1
    figure
    plot(unwrap(ang2(:,2)),'b-');
    %hold on;plot(unwrap(ang-phi_s),'r-');
    xlabel('step number')
    ylabel('unwrapped phase')
    figure
    scatter(rms,rmp,5,abs(ang2(:,1)));
    xlabel('rms')
    ylabel('rmp')
end
angfit = ang([30:70]);
%%
[ang_coeff, angp]=polyfit(xstep(fitrange),unwrap(ang2(fitrange,2)),1);
if plotflag == 1
    figure;
    plot(xstep(fitrange),unwrap(ang2(fitrange,2)),'b');
    %hold on
    %plot(xstep(fitrange),unwrap(ang(fitrange)-phi_s),'r');
    %plot(xstep(fitrange),xstep(fitrange)*ang_coeff(1)+ang_coeff(2),'--r');
    xlabel('relative stage position (nm)')
    ylabel('unwrapped phase')
end
end