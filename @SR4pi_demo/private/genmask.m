function [masksub] = genmask(res,lim)
maskLL = res.LL<lim.LL;
maskI = res.I>lim.I;
maskbg = res.bg<lim.bg;

maskz = res.z>lim.z(1) & res.z<lim.z(2);
masksubx = res.x>lim.x(1) & res.x<lim.x(2);
masksuby = res.y>lim.y(1) & res.y<lim.y(2);

masksub = masksubx & masksuby & maskLL & maskz & maskI & maskbg;

end
