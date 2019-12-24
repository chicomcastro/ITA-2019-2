function [eulerVetor] = quatToRotation(q)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% lambda = cos(phi/2) -> phi = 2*acos(lambda)
% rou_x = sin(phi/2)*phi_x/phi
q = q(:);

lambda = q(1);
r = q(2:4);

phi = 2*acos(lambda);
eulerVetor = r/sin(phi/2)*phi;

end

