function [deltaR] = CASO3_erroNavegacao(t,Omega,modTerra,deriva)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
e = deriva;
omega_s = modTerra.g_0/modTerra.R_0;

k = modTerra.g_0/((omega_s^2-Omega^2)*omega_s*Omega));

deltaR = [...
    0 k*(omega_s*sin(Ometa*t)-Omega*sin(omega_s*t)) -k/omega_s*(Omega^2*(1-cos(omega_s*t))-omega_s^2*(1-cos(Omega*t))); ...
    -modTerra.R_0/omega_s*(omega_s*t-sin(omega_s*t)) 0 0; ...
    0 0 0 ...
    ]*e;
end