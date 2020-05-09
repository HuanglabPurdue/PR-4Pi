function [shift1, shift2, shift3]=drift_correction_core3d_stack(im1,im2,iniguess)
%kernelsize=2;
if isempty(im1)||isempty(im2)||(sum(im1(:))==0)||(sum(im2(:))==0)
    shift1=1e25;
    shift2=1e25;
    shift3=1e25;
else
    corrim3d=normxcorr3_sparse(double(im1),double(im2),'same');
    % corrim3d=normxcorr3_sparse(double(gaussf(im1,kernelsize)),double(gaussf(im2,kernelsize)),'same');
    % corrim_norm=double(gaussf(normxcorr2(im1,im1),0));
    % scaleim=corrim-corrim_norm;
    % scaleim(~(scaleim<1e37&scaleim>-1e37))=0;
    
    if nargin>=3
        shift_col=iniguess(1);
        shift_row=iniguess(2);% this is the y axis and column in dipshow
        shift_z=iniguess(3);
        rowval=shift_row+round((size(corrim3d,1)+1)./2);
        colval=shift_col+round((size(corrim3d,2)+1)./2);
        zval=shift_z+round((size(corrim3d,3)+1)./2);
        cropsz2 = 24;
        rangelim_all = 0;
        while rangelim_all ~= 1
            ori = [colval-cropsz2/2-1 rowval-cropsz2/2-1 zval-cropsz2/2-1];
            smallim = cut(corrim3d,cropsz2,ori);
            %         data = double(smallim);
            %         P = zeros(5,size(data,3));
            %         for ii = 1:size(data,3)
            %             P(:,ii) = gaussfit(data(:,:,ii),[0.2,8,0.1,0.5,0.5]);
            %         end
            %
            %         x0 = -1.*P(4,:);
            %         y0 = -1.*P(5,:);
            %         I = 1;
            %         sigma = P(2,:);
            %         sigma([1:2,end-2:end]) = max(sigma);
            %         bg = 0;
            %         model = gengaussian(x0,y0,I,sigma,bg,cropsz2);
            %
            %         [tmpval colval2 rowval2 zval2] = findmax3d(double(model));
            mask = ones(size(smallim));
            mask(2:end-1,2:end-1,:) = 0;
            bg = zeros(1,1,size(smallim,3));
            bg(:,:,:) = squeeze(sum(sum(double(smallim.*mask),1),2))./sum(sum(mask(:,:,1)));
            bgL = repmat(bg,[cropsz2,cropsz2,1]);
            [tmpval colval2 rowval2 zval2] = findmax3d(double(smallim-bgL));
            colval=ori(1)+colval2;
            rowval=ori(2)+rowval2;
            zval=ori(3)+zval2;
            rangelim1 = colval2>cropsz2/4 && colval2<3*cropsz2/4;
            rangelim2 = rowval2>cropsz2/4 && rowval2<3*cropsz2/4;
            rangelim3 = zval2>cropsz2/4 && zval2<3*cropsz2/4;
            rangelim_all = rangelim1 && rangelim2 && rangelim3;
        end
    else
        [tmpval colval rowval zval] = findmax3d(double(corrim3d));
    end
    
    shift_row=(rowval-(size(corrim3d,1)+1)./2);
    shift_col=(colval-(size(corrim3d,2)+1)./2);
    shift_z=(zval-(size(corrim3d,3)+1)./2);
    cropsz=12;
    exsz=256;
    ori=[colval2-cropsz/2-1 rowval2-cropsz/2-1 zval2-cropsz/2-1];
    rscorrim=ift(extend(ft(cut(smallim-bgL,cropsz,ori)),exsz));
    [tmpval colval rowval zval]=findmax3d(double(abs(rscorrim)));
    s2row=(rowval-(exsz+1)/2)*cropsz/exsz;
    s2col=(colval-(exsz+1)/2)*cropsz/exsz;
    s2z=(zval-(exsz+1)/2)*cropsz/exsz;
    shift1=shift_col+s2col;
    shift2=shift_row+s2row;
    shift3=shift_z+s2z;
end