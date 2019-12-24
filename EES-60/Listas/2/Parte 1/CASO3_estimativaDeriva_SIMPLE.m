function [e] = CASO3_estimativaDeriva_SIMPLE(t,modTerra,erroHorizontal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R = raiosModTerra(modTerra,0);
R_e = R(3);
g_e = modTerra.g_0;

e = erroHorizontal/...
    sqrt(2*(R_e*t)^2 + 2*(R_e*sqrt(R_e/g_e)*sin(sqrt(g_e/R_e)*t))^2);

end