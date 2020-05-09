% for v14 iPALMast_analysisv14scmos.m

function smoothim=iPALM_unif_sCMOS(im,varim,sz)
if nargin<3
sz=3;
end
% imunf=unifim=varunif(im,varim,5)
smoothim=varunif(im,varim,sz)-varunif(im,varim,2*sz+3);
%-unif(im,[2*sz+3 2*sz+3 0],'rectangular');