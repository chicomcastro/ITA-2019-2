function [J] = indiceNormalidade(q_NED_B_comp)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    J = abs(norm(q_NED_B_comp)-1);

end

