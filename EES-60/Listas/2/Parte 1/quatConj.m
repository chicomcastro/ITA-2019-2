function [q_conj] = quatConj(q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

q = q(:);
lambda = q(1);
r = q(2:4);
q_conj = [lambda;-r(:)];

end

