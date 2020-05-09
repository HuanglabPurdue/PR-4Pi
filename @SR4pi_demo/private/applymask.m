function strout = applymask(strin,mask)

fd = fields(strin);
for ii = 1:numel(fd)
    strout.(fd{ii}) = strin.(fd{ii})(mask,:);
end

end