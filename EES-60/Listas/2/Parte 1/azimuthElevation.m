function [azimute,elevation] = azimuthElevation(lat_ref, lon_ref, h_ref, coordSV_e)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Constantes
Rt = 6378e3;
eixo1 = [1,0,0];
eixo2 = [0,1,0];
eixo3 = [0,0,1];

% Coord geodesicas (Terra esférica)
LAT_receptor = lat_ref;
LON_receptor = lon_ref;

% DCM
D_eNED = @(lat, lon) rotationToDCM(-(lat+pi/2)*eixo2)*rotationToDCM(lon*eixo3);

% ECEF (S_e)
X_ref_e = [...
    (Rt+h_ref)*cos(lat_ref)*cos(lon_ref); ...
    (Rt+h_ref)*cos(lat_ref)*sin(lon_ref); ...
    (Rt+h_ref)*sin(lat_ref) ...
    ];

% Azimute e elevação
r_e = coordSV_e(:) - X_ref_e(:);  % em S_e
r_NED = D_eNED(LAT_receptor,LON_receptor)*r_e/norm(r_e);  % norm no NED do receptor

rn = r_NED(1);
re = r_NED(2);
rd = r_NED(3);

azimute = rad2deg(mod(atan2(re,rn),2*pi));
elevation = rad2deg(atan(-rd/sqrt(rn^2+re^2)));

end
