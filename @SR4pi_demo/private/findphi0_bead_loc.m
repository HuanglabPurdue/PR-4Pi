
% find phic, phase shift at z=0
function [phic,phic_u] = findphi0_bead_loc(sr,phi,obj)
srcenter = obj.Srcenter;
phi_offset = obj.Phi_offset;
plotflag = 0;
phia = unwrap(phi);
masksub = sr>-2&sr<2;
p_phi = polyfit(sr(masksub),phia(masksub),1);
phic = polyval(p_phi,srcenter);
if plotflag ==1
    figure;
    plot(sr(masksub),phia(masksub),'.')
    hold on
    plot(srcenter,phic,'o')
end
phic_u = phic;
phic = wrapToPi(phic+phi_offset);
end
