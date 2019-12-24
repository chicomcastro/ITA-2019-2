%%
% Coordenadas dadas no ECEF [m]
XYZ = [...
    7766188.44, -21960535.34, 12522838.56; ...  % SV2
    -25922679.66, -6629461.28, 31864.37; ...    % SV26
    -5743774.02, -25828319.92, 1692757.72; ...  % SV4
    -2786005.69, -15900725.80, 21302003.49 ...  % SV7
    ];

X(:) = XYZ(:,1); X = X(:);
Y(:) = XYZ(:,2); Y = Y(:);
Z(:) = XYZ(:,3); Z = Z(:);
%%
% Pseudodistâncias [m]
rho(1) = 22228206.42;
rho(2) = 24096139.11;
rho(3) = 21729070.63;
rho(4) = 21259581.09;
rho = rho(:);

%A distância computada ao i-ésimo satélite é i i
%c deltaT ? c ?tau deltar sendo 

%%
% Estimativa inicial
X0 = 0;
Y0 = 0;
Z0 = 0;
b0 = 0;

deltaX = 0;

i = 0;
iteracoes(1,:) = [i X0 0 Y0 0 Z0 0 b0];
while (1)

    rho_chapeu = ((X-X0).^2+(Y-Y0).^2+(Z-Z0).^2).^0.5 + b0;

    a(:,1) = (X0 - X)./(rho_chapeu - b0);
    a(:,2) = (Y0 - Y)./(rho_chapeu - b0);
    a(:,3) = (Z0 - Z)./(rho_chapeu - b0);
    a(:,4) = 1;

    deltarho = rho - rho_chapeu;

    deltaX = inv(a)*deltarho;
    X0 = X0 + deltaX(1);
    Y0 = Y0 + deltaX(2);
    Z0 = Z0 + deltaX(3);
    b0 = b0 + deltaX(4);
    erro = deltaX;
    
    iteracoes(i+1,:) = [i+1 X0 deltaX(1) Y0 deltaX(2) Z0 deltaX(3) b0];
    i = i+1;
    
    if(norm(erro) < 1e-8)
        break;
    end
end

X_e = [X0;Y0;Z0;b0];

%%
% Constantes
eixo1 = [1,0,0];
eixo2 = [0,1,0];
eixo3 = [0,0,1];

% Coord geodesicas (Terra esférica)
LAT_receptor = asin(Z0/norm([X0 Y0 Z0]));
LON_receptor = atan2(Y0,X0);

% DCM
D_eNED = @(lat, lon) rotationToDCM(-(lat+pi/2)*eixo2)*rotationToDCM(lon*eixo3);

% Azimute e elevação
r_e = XYZ' - X_e(1:3);  % em S_e
r_NED = D_eNED(LAT_receptor,LON_receptor)*r_e/norm(r_e);  % norm no NED do receptor
for i = 1:4
    rn = r_NED(1,i);
    re = r_NED(2,i);
    rd = r_NED(3,i);
    
    azimute(i) = rad2deg(mod(atan2(re,rn),2*pi));
    elevation(i) = rad2deg(atan(-rd/sqrt(rn^2+re^2)));
end


%% Plot

axis = polaraxes;
for i = 1:4
    polarplot(axis, deg2rad(azimute(i)), elevation(i), 'h','LineWidth',2);
    hold on;
end

axis.ThetaDir = 'clockwise';
axis.ThetaZeroLocation = 'top';

axis.RLim = [0 90];
axis.RDir = 'reverse';
grid on

legend('SV 2','SV 26','SV 4','SV 7');
