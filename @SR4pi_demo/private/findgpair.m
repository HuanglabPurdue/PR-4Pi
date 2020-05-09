function [gpair,gzone] = findgpair(res_fitsub,gzone,plotflag)
disx = pdist(cat(2,res_fitsub.x,res_fitsub.y),'euclidean');
disz = pdist(cat(2,res_fitsub.z,zeros(size(res_fitsub.z))),'euclidean');
if isempty(gzone.xlim)
    figure;
    ha = axes;
    plot(disx,disz,'.','markersize',1);axis equal;
    ha.XLim = [0,100];
    ha.YLim = [200,400];
    drawnow;
    rect1 = getrect(ha);
    gzone.xlim = [rect1(1),rect1(1)+rect1(3)];
    gzone.zlim = [rect1(2),rect1(2)+rect1(4)];
    gzone.zT = mean(gzone.zlim);
end
N = numel(res_fitsub.x);
maskg = int8(disx<gzone.xlim(2) & disz>gzone.zlim(1) & disz<gzone.zlim(2));
ind = find(maskg==1);
Ni = numel(ind);
indL = repmat(ind',1,N-1);
a = [1:N-1];
lowM = repmat((a-1).*N-(a-1).*a./2+1,Ni,1);
highM = repmat((a-1).*N-(a-1).*a./2+N-a,Ni,1);

tmp1 = indL-lowM;
tmp2 = highM-indL;
maski = int8(tmp1>=0 & tmp2>=0);
[~,inda] = find(maski==1);
indb = ind'-(inda-1).*N+(inda-1).*inda./2+inda;


tmp = zeros(N,1);
tmp(inda) = 1;
indra = find(tmp==1);

Nt = numel(indra);
gpair = cell(Nt,1);

for ii = 1:Nt
    T = inda == indra(ii);
    gpair{ii} = cat(1,indra(ii),indb(T));
end

for nn = 1:3
    for kk = 1:10
        count = 1;
        for ii = 1:Nt
            if size(gpair,1)>count+kk-1
                tmp = cat(1,gpair{count},gpair{count+kk});
                tmp1 = unique(tmp);
                if numel(tmp)>numel(tmp1)
                    gpair{count} = tmp1;
                    gpair(count+kk) = [];
                else
                    count = count + 1;
                end
            else
                break;
            end
        end
    end
end

if plotflag == 1
    figure;
    ha = axes;
    hold(ha,'on');
    for ii = 1:size(gpair,1)
        plot(res_fitsub.y(gpair{ii}),res_fitsub.z(gpair{ii}),'o','parent',ha)
    end
    axis equal
end
