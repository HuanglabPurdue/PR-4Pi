% for quadrant alignment

function [q1, q2, q3, q4] = align(qd1,qd2,qd3,qd4,affS)
zm = affS.zm_all;
trans = affS.trans_all;
ang = affS.ang_all;

parfor ii = 1:size(qd1,3)
    
    q1(:,:,ii) = qd1(:,:,ii);
    
    q2(:,:,ii) = double(affine_trans(qd2(:,:,ii),zm(2,:),trans(2,:),ang(2)));
    
    q3(:,:,ii) = double(affine_trans(qd3(:,:,ii),zm(3,:),trans(3,:),ang(3)));
    
    q4(:,:,ii) = double(affine_trans(qd4(:,:,ii),zm(4,:),trans(4,:),ang(4)));
    
end


