function [q] = rotationToQuat(phi,dir)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

q = [
    cos(phi/2);
    sin(phi/2)*sin(dir(1)*pi/2);
    sin(phi/2)*sin(dir(2)*pi/2);
    sin(phi/2)*sin(dir(3)*pi/2)
    ];

end

