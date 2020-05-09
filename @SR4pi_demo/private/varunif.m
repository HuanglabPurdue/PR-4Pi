
function unifim=varunif(im,varmap,sz)
im=single(im);
varmap=single(varmap);
if size(varmap,3)==1
    varmap=repmat(varmap,[1 1 size(im,3)]);
end
tmp1=im./varmap;
kerim=ones(sz,sz);
cim=convn(tmp1,kerim,'same');
wim=convn(1./varmap(:,:,1),kerim,'same');
wim=repmat(wim,[1 1 size(cim,3)]);
unifim=cim./wim;