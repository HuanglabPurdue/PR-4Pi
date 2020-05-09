function out=mephigrow(mephi,dmap,maxangle,srange,stopval,invflag,ha,ha1)
I=1e37;
msz=size(dmap,1);
while I>stopval
    v1=(mephi(end,1)-mephi(end-1,1))+wrapToPi(mephi(end,2)-mephi(end-1,2))*1i;
    %     v2=(coords(3,1)-coords(2,1))+(coords(3,2)-coords(2,2))*1i;
    %     phi=wrapToPi(acos(((coords(2,1)-coords(1,1))*(coords(3,1)-coords(2,1))+(coords(2,2)-coords(1,2))*(coords(3,2)-coords(2,2)))./abs(v1)./abs(v2)));
    xx=linspace(-pi,pi,msz);
    yy=linspace(-pi,pi,msz);
    [xgrid ygrid]=meshgrid(xx,yy);
    rmap=sqrt((xgrid-mephi(end,1)).^2+wrapToPi(ygrid-mephi(end,2)).^2);
    anglemap=wrapToPi(acos(((mephi(end,1)-mephi(end-1,1)).*(xgrid-mephi(end,1))+...
        wrapToPi(mephi(end,2)-mephi(end-1,2))*wrapToPi(ygrid-mephi(end,2)))./abs(v1)./abs(rmap)));% inner product of two vector, angel=0 or pi along the vector formed by the last two points of mephi
    
    maska=anglemap<maxangle;
    maskr=rmap<srange(2)&rmap>srange(1);
    directmask=(xgrid-mephi(end,1))./(mephi(end,1)-mephi(end-1,1))>=invflag.*0.5;
    
    mask=maska.*maskr.*directmask;
    
    vismap=mask.*dmap;
    tmpmap = ha.UserData + vismap;
    imagesc(tmpmap,'parent',ha);
    ha.UserData = tmpmap;
%     keyboard
    [I row col]=findmax(vismap);
    mephi(end+1,2)=yy(row);
    mephi(end,1)=xx(col);
    
    scatter(ha1,mephi(:,1),mephi(:,2));
end

out=mephi;