function [rch gch bch]=srhist_color(sz,zm,xtot,ytot,ttot,segnum,colorstring)
% [rch gch bch]=srhist_color(sz,zm,xtot,ytot,ttot,segnum,colorstring)
if nargin<7
    colorstring='hsv';
end
ttot=ttot(1:size(xtot,1));
eval(['a=colormap(' colorstring '(segnum+20));']);
a = a(21:end,:);
%close all
% r = [255,253,252,250,247,221,174,122,73];
% g = [247,224,197,159,104,52,1,1,0];
% b = [243,221,192,181,161,151,126,119,106];

% r = [215   244   253   254   255   224   171   116    69];
% g = [48   109   174   224   255   243   217   173   117];
% b = [39    67    97   144   191   248   233   209   180];

% r = [197   222   241   253   247   230   184   127    77];
% g = [27   119   182   224   247   245   225   188   146];
% b = [125   174   218   239   247   208   134    65    33];

% r = [255   232    93     0     0     0     0    93   232];
% g = [139   255   255   255   255   185    46     0     0];
% b = [0     0     0    46   185   255   255   255   255];
% 
% rL = interp1([0:8:64],r,[0:64]);
% gL = interp1([0:8:64],g,[0:64]);
% bL = interp1([0:8:64],b,[0:64]);
% 
% a = flip(cat(1,rL,gL,bL)'./255);
load('cmap6.mat')
a = cm;
incre=floor((max(ttot)-min(ttot)+1)/segnum);
for ii=1:1:segnum
    tst=(ii-1)*incre+min(ttot);
    if ii==segnum
        ted=max(ttot);
    else
        ted=ii*incre++min(ttot);
    end
    
    mask=ttot>=tst&ttot<=ted;
    xtmp=xtot(mask);
    ytmp=ytot(mask);
    
    tmpim=SRreconstructhist(sz,zm,xtmp,ytmp);
    tmpim=double(tmpim);
    if ii==1
        rch=a(ii,1)*tmpim;
        gch=a(ii,2)*tmpim;
        bch=a(ii,3)*tmpim;
    else
        rch=rch+a(ii,1)*tmpim;
        gch=gch+a(ii,2)*tmpim;
        bch=bch+a(ii,3)*tmpim;
    end
end
