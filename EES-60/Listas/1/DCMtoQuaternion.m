function [q] = DCMtoQuaternion(DCM)

q = [
    sqrt((trace(DCM)+1)/4);
    (DCM(2,3)-DCM(3,2))/(4*sqrt((trace(DCM)+1)/4));
    (-DCM(1,3)+DCM(3,1))/(4*sqrt((trace(DCM)+1)/4));
    (DCM(1,2)-DCM(2,1))/(4*sqrt((trace(DCM)+1)/4))
    ];

end

