function [R_sat_e,V_sat_e_inercial] = coordSV(SV,t,Sat_OrbitData)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Constantes
eixo1 = [1,0,0];
eixo2 = [0,1,0];
eixo3 = [0,0,1];

% Dados
Omega = 7.29e-5;    % rad/s
R = 26560000;       % [m]
T = 43.082;         % [s]
i = 55*pi/180;      % rad
mi = 3.986e14;

% o argumento de latitude u(t) varia com taxa constante
udot = 2*pi/T;
u = @(t, u0) u0 + udot*t;

% Em vez da ascensão reta do nó ascendente (RAAN), será usada a longitude geográfica do nó
% ascendente glan(t), cuja taxa de variação constante se deve à velocidade angular ? da
% Terra
glan = @(t, glan0) glan0 - Omega*t;

D_e_O = @(t, u0, glan0) rotationToDCM( u(t, u0)*eixo3 )*...
    rotationToDCM( i*eixo1 )*...
    rotationToDCM( glan(t, glan0)*eixo3 );

%%
% Posição e velocidade inercial dos 24 SVs no ECEF
R_sat_O = [R 0 0]';
V_sat_O = [0 sqrt(mi/R) 0]';

% Dado um satélite (SV) e um tempo, R_sat_e e V_sat_e_inercial são

u0 = Sat_OrbitData(SV,1);
glan0 = Sat_OrbitData(SV,2);
D_O_e = inv(D_e_O(t,u0,glan0)); % D_O_e
R_sat_e = D_O_e*R_sat_O;
V_sat_e = D_O_e*V_sat_O;

Omega_e = [0,0,Omega];
V_sat_e_inercial = V_sat_e(:) + cross(Omega_e(:),R_sat_e(:));

end

