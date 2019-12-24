function [q_inv] = quatInv(q)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
q = q(:);
q_inv = quatConj(q)/norm(q);

end

