function [subroi,centers]=pickbead(data,boxsize,N,type,centers)
s1 = floor(boxsize/2);
s2 = ceil(boxsize/2);

sz = size(data,3);
xsz = max(size(data,1),size(data,2));
if isempty(centers)
    h=dipshow(data(:,:,round(sz/2)));
    diptruesize(h,round(1000/xsz*100))
    centers = dipgetcoords(h,N);
    close(h)
end

switch type
    case 'single'
        subroi=zeros(boxsize,boxsize,N);
        for ii=1:N
            roi=data(centers(ii,2)-s1+1:centers(ii,2)+s2,centers(ii,1)-s1+1:centers(ii,1)+s2,centers(ii,3)+1);
            subroi(:,:,ii)=roi;
        end
    case 'stack'
        subroi=[];
        for ii=1:N
            roi=data(centers(ii,2)-s1+1:centers(ii,2)+s2,centers(ii,1)-s1+1:centers(ii,1)+s2,:,:);
            subroi=cat(3,subroi,roi);
        end
end

