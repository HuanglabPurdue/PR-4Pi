% for v14 iPALMast_analysisv14scmos.m
% Use astfun_val to give astigmatism back

function [estimates sse]=iPALMast_iniz_astfit(sx,sy,astparam)
ini=1000;
model = @Ast_est;
[estimates sse] = fminsearch(model, ini);
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = Ast_est(params)
        z0 = params(1);
        exsigx = astfun_val(z0,astparam.estx);
        exsigy = astfun_val(z0,astparam.esty);
        ErrorVector = sqrt((exsigx-sx).^2+(exsigy-sy).^2)+1000.*abs(imag(exsigx))+1000.*abs(imag(exsigy));
        sse = sum(ErrorVector(:));
    end
end