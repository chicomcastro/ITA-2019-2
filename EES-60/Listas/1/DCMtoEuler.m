function [euler] = DCMtoEuler(DCM)
%UNTITLED3 Summary of this function goes here
euler = [
    atan(DCM(1,2)/DCM(1,1));
    atan(-DCM(1,3)/(1-DCM(1,3)^2)^(0.5));
    atan(DCM(2,3)/DCM(3,3))
    ];

end

