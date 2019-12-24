function [DCM] = quatToDCM(q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
q = q(:);
DCM = [
    q(1)^2+q(2)^2-q(3)^2-q(4)^2     2*(q(2)*q(3)+q(1)*q(4))         2*(q(2)*q(4)-q(1)*q(3));
    2*(q(2)*q(3)-q(1)*q(4))         q(1)^2-q(2)^2+q(3)^2-q(4)^2     2*(q(3)*q(4)+q(1)*q(2));
    2*(q(2)*q(4)+q(1)*q(3))         2*(q(3)*q(4)-q(1)*q(2))         q(1)^2-q(2)^2-q(3)^2+q(4)^2
    ];

end

