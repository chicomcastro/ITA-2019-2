%% Lab 1 - EES-60 - versão 2

%% Valores envolvidos
tetha_c = 11.4/180*pi; % rad
Omega_prec = 30/180*pi; % rad/s

%% Solução analítica

psi = @(t) Omega_prec*t;
tetha = @(t) tetha_c;
phi = @(t)  -Omega_prec*t;


%% Quaternions de rotação

q1 = @(t) rotationToQuat(psi(t),[0,0,1]);
q2 = @(t) rotationToQuat(tetha(t),[0,1,0]);
q3 = @(t) rotationToQuat(phi(t),[0,0,1]);

q_NED_B = @(t) quatProd(quatProd(q3(t),q2(t)),q1(t));


%% Velocidade angular do corpo em relação ao NED descrito na base do corpo
w_B_NED_x = @(t) -Omega_prec*sin(tetha_c)*cos(Omega_prec*t);
w_B_NED_y = @(t) -Omega_prec*sin(tetha_c)*sin(Omega_prec*t);
w_B_NED_z = @(t) Omega_prec*(cos(tetha_c)-1);


%% Matriz de relação dos quaternions para a EDO
Omega = @(t)[
    0             -w_B_NED_x(t) -w_B_NED_y(t) -w_B_NED_z(t);
    w_B_NED_x(t)  0             w_B_NED_z(t)  -w_B_NED_y(t);
    w_B_NED_y(t)  -w_B_NED_z(t) 0             w_B_NED_x(t);
    w_B_NED_z(t)  w_B_NED_y(t)  -w_B_NED_x(t) 0
    ];


%% EDO
f = @(t,y) 1/2*Omega(t)*y;              % y' = f(t,y)
y = q_NED_B(0);          % y(t0) = y0
h = 1e-1;                               % passo

% Condições de parada
tf = 100;                               % tempo final (s)


%% Data handlers
Psi = zeros(tf/h,3);
result = zeros(tf/h,4);


%% Integração numérica

i = 1;

for t=0:h:tf
    %% Rughe-kutta passo fixo 4a ordem
    k1 = f(t,y);
    k2 = f(t + h/2, y + h/2*k1);
    k3 = f(t + h/2, y + h/2*k2);
    k4 = f(t + h, y + h*k3);

    y = y + h/6*(k1 + 2*k2 + 2*k3 + k4);

    result(i,:) = y;
    
    %% Desalinhamento
    y = y/norm(y);
    D_NED_B_comp = quatToDCM(y);
    D_NED_B_anal = quatToDCM(q_NED_B(t));
    D_NED_NED_comp = D_NED_B_comp'*D_NED_B_anal;

    %% Vetor desalinhamento (em radianos)
    Psi(i,:) = [
        ( abs(D_NED_NED_comp(2,3)) + abs(D_NED_NED_comp(3,2)) ) / 2;
        ( abs(D_NED_NED_comp(1,3)) + abs(D_NED_NED_comp(3,1)) ) / 2;
        ( abs(D_NED_NED_comp(1,2)) + abs(D_NED_NED_comp(2,1)) ) / 2
        ];
    i = i+1;
end

%% Plots
figure;
plot(0:h:tf, Psi); % Plotar em arcseg <<<<<

% O que é o ângulo de desalinhamento em torno do eixo de desalinhamento?
% <<<<

%% Verificação
figure;
plot(0:h:tf, result);
hold on;
[TOUT,YOUT] = ode45(f, [0 tf], q_NED_B(0));
plot(TOUT,YOUT,'.');
    
