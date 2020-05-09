% t: vector of frame numbers
% nfs: number of frame per section

function [cormask] = gencormask(t,nfs)
fn = ceil(single(max(t))/nfs);
v = [0:fn];
cormask = zeros(size(t));
for ii = 1:numel(v)-1
    mask = (t >= v(ii)*nfs)&(t < v(ii+1)*nfs);
    cormask(mask) = ii;
end
cormask(cormask==0) = ii;
end