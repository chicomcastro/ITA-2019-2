
A = 4000; % kg.m^2
B = 7500; % kg.m^2
C = 8500; % kg.m^2

I = [
    A 0 0;
    0 B 0;
    0 0 C
    ];
M = [0,0,0]; % N.m

y0 = [...
    0.1;...
    -0.2;...
    0.5;...
    0;...
    pi\2;...
    0 ...
    ];

[t,y] = ode45(@(t,y) f(t,y,I,M), [0 100], y0);

resultados.tempo = t;
resultados.omega = y(:,1:3);
resultados.atitude = y(:,4:6);

figure;
plot(resultados.tempo, resultados.omega)
title("Velocidade angular");
ylabel('Velocidade angular [rad/s]');
xlabel('Tempo [s]');
legend('\omega_x', '\omega_y', '\omega_z');

figure;
plot(resultados.tempo, resultados.atitude);
ylabel('Atitude [rad]');
xlabel('Tempo [s]');
legend('yaw', 'pitch', 'roll');


function [y_ponto] = f(t,y,I,M)

A = I(1,1);
B = I(2,2);
C = I(3,3);
M = M(:);
M_x = M(1);
M_y = M(2);
M_z = M(3);
omega_x = y(1);
omega_y = y(2);
omega_z = y(3);
phi = y(4);
theta = y(5);
psi = y(6);

y_ponto = [...
    (M_x-(C-B)*omega_y*omega_z)/A;...
    (M_y-(A-C)*omega_z*omega_x)/B;...
    (M_z-(B-A)*omega_x*omega_y)/C;...
    (omega_x*sin(psi)+omega_y*cos(psi))/sin(theta);...
    omega_x*cos(psi)-omega_y*sin(psi);...
    -(omega_x*sin(psi)+omega_y*cos(psi))/tan(theta)+omega_z...
    ];

end