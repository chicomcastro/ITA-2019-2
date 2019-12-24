%% Questão 03

% Abaixo se encontram, respectivamente, os valores
% iniciais do argumento de latitude u(0) e da longitude geográfica do nó ascendente glan(0)
% na época 01/julho/1990 à 00h00’00” (Table I. Orbit Reference Values of the 1993 edition
% of the Navstar GPS description report by NavtechGPS).
Sat_OrbitData(1:24,1:2)=0;  % initial argument of latitude u(0) and
                            % geographic longitude of ascending node
                            % glan(0) [rad]
Sat_OrbitData(1,1:2)=[280.7 358.6]*pi/180;     % SVs 1 to 10
Sat_OrbitData(2,1:2)=[310.3 13.4]*pi/180;
Sat_OrbitData(3,1:2)=[60 248.2]*pi/180;
Sat_OrbitData(4,1:2)=[173.4 304.9]*pi/180;
Sat_OrbitData(5,1:2)=[339.7 88.1]*pi/180;
Sat_OrbitData(6,1:2)=[81.9 319.2]*pi/180;
Sat_OrbitData(7,1:2)=[115.0 335.7]*pi/180;
Sat_OrbitData(8,1:2)=[213.9 25.2]*pi/180;
Sat_OrbitData(9,1:2)=[16 346.2]*pi/180;
Sat_OrbitData(10,1:2)=[138.7 47.6]*pi/180;
Sat_OrbitData(11,1:2)=[244.9 100.7]*pi/180;    % SVs 11 to 20

Sat_OrbitData(12,1:2)=[273.5 115]*pi/180;
Sat_OrbitData(13,1:2)=[42.1 59.3]*pi/180;
Sat_OrbitData(14,1:2)=[70.7 73.6]*pi/180;
Sat_OrbitData(15,1:2)=[176.8 126.6]*pi/180;
Sat_OrbitData(16,1:2)=[299.6 188]*pi/180;
Sat_OrbitData(17,1:2)=[101.7 149.1]*pi/180;
Sat_OrbitData(18,1:2)=[200.5 198.5]*pi/180;
Sat_OrbitData(19,1:2)=[233.7 215.1]*pi/180;
Sat_OrbitData(20,1:2)=[335.9 266.2]*pi/180;

Sat_OrbitData(21,1:2)=[142.2 229.3]*pi/180;    % SVs 21 to 24
Sat_OrbitData(22,1:2)=[255.6 286]*pi/180;
Sat_OrbitData(23,1:2)=[5.3 160.9]*pi/180;
Sat_OrbitData(24,1:2)=[34.8 175.6]*pi/180;



%%
% Dada latitude, longitude, altitude e o intervalo de tempo escolhido
% --> computar quais SVs estão visíveis;
% --> a evolução da velocidade relativa dos SVs visíveis em relação ao
% receptor ao longo das respectivas LOSs em [m/s]; e
% --> a magnitude da velocidade angular de cada LOS para SV que esteja 
% visível em [?/h];

eixo1 = [1,0,0];
eixo2 = [0,1,0];
eixo3 = [0,0,1];
Rt = 6378e3;
Omega = 7.29e-5;

% Geográfica
lat = 0;        % [rad]
lon = 0;        % [rad]
h = 600;        % [m]
tf = 20;
t0 = 0;

X_ref_e = [...
    (Rt+h)*cos(lat)*cos(lon); ...
    (Rt+h)*cos(lat)*sin(lon); ...
    (Rt+h)*sin(lat) ...
    ];
            
D_NED_LOS = @(elevation, azimuth) rotationToDCM(elevation*eixo2)*rotationToDCM(azimuth*eixo3);
D_e_NED = @(lat, lon) rotationToDCM(-(lat+pi/2)*eixo2)*rotationToDCM(lon*eixo3);

m = (1+3+3+1+1+1+2);
result = zeros(m*24,tf-t0);
index = 1;

for t = t0:tf
    for i = 1:24
        [isVisible, R_sat_e, V_sat_e_inercial, azimuth, elevation] = isSVVisibleFrom([lat;lon;h],i,t,Sat_OrbitData);
        
        r_rel_e = R_sat_e(:) - X_ref_e(:);
        v_rel_e = V_sat_e_inercial(:) - cross([0;0;Omega],X_ref_e(:));

        v_rel_los = D_NED_LOS(deg2rad(elevation),deg2rad(azimuth))*D_e_NED(lat,lon)*v_rel_e;

        result(m*(i-1)+1:m*(i-1)+m, index) = [...
            t; ...
            R_sat_e(:); ...
            V_sat_e_inercial(:); ...
            isVisible; ...
            v_rel_los(1); ...
            norm( [0;v_rel_los(2:3)]/norm(r_rel_e) )*180/pi*3600; ...
            azimuth; ...
            elevation ...
        ];
    end
    index = index + 1;
end

%%
% Fazer skyplot da evolução dos azimutes e elevações dos SVs visíveis no
% referido intervalo
azimutes = [];
elevacoes = [];
visibilidades = [];
legenda = string();

figure;

axis = polaraxes;

for i = 1:24
    visibilidades = [visibilidades; result(m*(i-1)+8,:)];
    azimutes = [azimutes; result(m*(i-1)+11,:)];
    elevacoes = [elevacoes; result(m*(i-1)+12,:)];

    polarplot(axis, deg2rad(result(m*(i-1)+11,:)), result(m*(i-1)+12,:),'-');

    if(any(result(m*(i-1)+8,:)>0))
    legenda(end+1) = "SV " + i;
    end
    hold on
end

axis.ThetaDir = 'clockwise';
axis.ThetaZeroLocation = 'top';

axis.RLim = [0 90];
axis.RDir = 'reverse';
grid on

legend(legenda(2:end));
