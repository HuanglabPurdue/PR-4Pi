
function [est, sse] = iniz_astfit(sx,sy,astparam,z0)
if nargin>3
    ini = z0;
else
    zs = [-1,-0.5,0.5,1];
    %zs = 0.5*ones(1,4);
    sse0 = zeros(1,4);
    for ii = 1:4
        [sse0(ii)] = ast_est(zs(ii),astparam,sx,sy);
    end
    [~,ind] = min(sse0);
    ini = zs(ind);
end
%ini = 1;
est = fminsearch(@(x)ast_est(x,astparam,sx,sy),ini);%,optimset('MaxIter',1000,'Display','off','TolX',1e-12,'TolFun',1e-12));
% opts = optimoptions('fminunc','MaxIter',100,'algorithm','quasi-newton');
% est = fminunc(@(x)ast_est(x,astparam,sx,sy),ini,opts);

[sse] = ast_est(est,astparam,sx,sy);    
end

function [sse] = ast_est(z0,astparam,sx,sy)
[~,exsigx,~] = astfuneval(astparam.estx,z0);
[~,~,exsigy] = astfuneval(astparam.esty,z0);
%ErrorVector = sqrt((exsigx-sx).^2+(exsigy-sy).^2);
ErrorVector = sqrt((exsigx-sx).^2+(exsigy-sy).^2)+1000.*abs(imag(exsigx))+1000.*abs(imag(exsigy));
%ErrorVector = ((exsigx-sx).^2+(exsigy-sy).^2)+1000.*abs(imag(exsigx))+1000.*abs(imag(exsigy));
sse = double(sum(ErrorVector(:)));
end