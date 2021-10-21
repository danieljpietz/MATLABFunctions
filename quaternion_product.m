function [qab] = quaternion_product(qa,qb)
qab = [qa(1) -qa(2) -qa(3) -qa(4);
    qa(2) qa(1) -qa(4) qa(3);
    qa(3) qa(4) qa(1) -qa(2);
    qa(4) -qa(3) qa(2) qa(1)] * qb;
end

