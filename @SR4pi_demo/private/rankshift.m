function [rshift,errorf]= rankshift(AA,shift,error,thresh)
newAA = AA;
oldAA = newAA;
oldshift = shift;
while rank(newAA)>=(min(size(AA))-1)&&max(error)>thresh
    oldAA = newAA;
    oldshift = shift;
    [~, ind]=max(error);
    newAA(ind,:) = [];
    error(ind) = [];
    shift(ind) = [];
end

if rank(newAA)<(min(size(AA))-1)
    newAA = oldAA;
    shift = oldshift;
    flag(1)=1;
end

rshift = pinv(newAA)*shift;
errorf = sqrt(((newAA*rshift)-shift).^2);
