function [deltaTheta] = erroPosVerdadeira(deltaR, coord, modTerra)
% Retorna deltaTheta em arcseg
N = 1;
E = 2;

deltaR = deltaR(:);

h = coord.alt;
lambda = coord.lat;
[R_N, R_E] = raiosModTerra(modTerra,lambda);

deltaTheta = [...
    deltaR(E)/(R_E+h); ...
    -deltaR(N)/(R_N+h); ...
    -deltaR(E)/(R_E+h)*tan(lambda) ...
    ]*180/pi*3600;
end

