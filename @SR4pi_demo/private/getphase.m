function [ang] = getphase(P,subims,plotflag)
%% 
xco=P(:,2);
yco=P(:,1);

%% get rms rmp
xf=xco;
yf=yco;
[rms, rmp]=iPALMast_findmom_givenC(subims,xf,yf);

%% calculate reduced moment
rm1=rms./sqrt(rms.^2+rmp.^2);
rm2=rmp./sqrt(rms.^2+rmp.^2);
c=rm1-1i*rm2;
ang=angle(c);   % change to MLE
if plotflag == 1
    figure;
    plot(ang);
end


end