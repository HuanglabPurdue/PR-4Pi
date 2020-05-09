function [ndata,mu,sigma] = normalizedist(data,value,plotflag)
if isscalar(value)
    [count,binc] = hist(data,value);
else
    binc = value;
    count = hist(data,binc);
end
f = fit(binc(2:end-1)',count(2:end-1)','gauss1');
if plotflag == 1
    figure;
    plot(f,binc(2:end-1),count(2:end-1))
end
mu = f.b1;
sigma = f.c1/sqrt(2);
ndata = (data - f.b1).*sqrt(2)./f.c1;

end