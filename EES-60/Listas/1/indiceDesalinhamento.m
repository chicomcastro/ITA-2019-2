function [I] = indiceDesalinhamento(A,B)
%UNTITLED5 Summary of this function goes here
%   B is the reference D0

    A = quatToDCM(A);
    B = quatToDCM(B);
    I = abs(acos(1/2*trace(A'*B)-1/2));

end

