function [render_img] = renderslice1(resi,xsz,zm,pixelsize,pmax,sigma,cmp)
% x,y in pixels, z in nm
[rch, gch, bch]=srhist_color(xsz,zm,resi.x./pixelsize,resi.y./pixelsize,resi.z,64,cmp);

low_out_z=0;
high_out_z=255;

rchsm = gaussf(rch,sigma);
gchsm = gaussf(gch,sigma);
bchsm = gaussf(bch,sigma);
normf = min([max(rchsm),max(gchsm),max(bchsm),pmax]);
ratio_r = min([max(rchsm),pmax])/normf;
ratio_g = min([max(gchsm),pmax])/normf;
ratio_b = min([max(bchsm),pmax])/normf;

rchsmst = imstretch_linear(rchsm,0,pmax*ratio_r,low_out_z,high_out_z*ratio_r);
gchsmst = imstretch_linear(gchsm,0,pmax*ratio_g,low_out_z,high_out_z*ratio_g);
bchsmst = imstretch_linear(bchsm,0,pmax*ratio_b,low_out_z,high_out_z*ratio_b);
colorim = joinchannels('RGB',rchsmst,gchsmst,bchsmst);

render_img = permute(double(colorim),[2,1,3]);
end

