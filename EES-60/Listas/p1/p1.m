eixo1 = [1,0,0];
eixo2 = [0,1,0];
eixo3 = [0,0,1];


%% i
syms alpha i u positive
D_IO = rotationToDCM(u*eixo3)*rotationToDCM(i*eixo1)*rotationToDCM(alpha*eixo3);

%% ii
syms OmegaT positive
D_Ie = rotationToDCM(OmegaT*eixo3);

%% iii
syms tri lambda positive
D_eNED = rotationToDCM(-(lambda+pi/2)*eixo2)*rotationToDCM(tri*eixo3);

%% iv
syms x gama positive
D_NEDLOS = rotationToDCM(gama*eixo2)*rotationToDCM(x*eixo3);

%% Questão 02
% satelite
syms a mi positive
R_sat_O = [a 0 0]';
V_sat_O = [0 sqrt(mi/a) 0]';
% veleiro
syms h R positive
R_vlr_NED = [0 0 -(R+h)]';
V_vlr_NED = [0 0 0]';

R_sat_NED = D_eNED * D_Ie * D_IO' * R_sat_O;
%%
r = R_sat_NED - R_vlr_NED;
r = r/norm(r);
azimute = acos(dot((r-dot(r,[0,0,1]')),[1,0,0]'));
elevation = acos(dot((r-dot(r,[0,0,1]')),r));

%% Questão 06

mi = 3.986e14;
R = 6378e3;
Omega = 7.292115e-5;
a = 26500e3;

T = 3600;
alpha = 60*pi/180;
i = 55*pi/180;
u = 200*pi/180;
lambda = -23*pi/180;
tri = -43*pi/180;
h = 0;
OmegaT = Omega*T;

%%
azimute = double(subs(azimute))*180/pi
elevation = double(subs(elevation))*180/pi

%% Questão 04

omega_O = D_IO * D_Ie * D_eNED * [Omega*cos(lambda),0,-Omega*sin(lambda)]';
V_sat_O_inercial = V_sat_O;
V_sat_O = V_sat_O_inercial - cross(omega_O,R_sat_O);

V_sat_NED = D_eNED * D_Ie * D_IO' * V_sat_O;

double(subs(V_sat_NED*dot(V_sat_NED/norm(V_sat_NED),r)))
