function trueCovAnalysis(T, k, multiFactor)

% Modelo verdadeiro
% x_ponto = F x + w
% w ~ norm(0,Q)
% y = H x + v
% v ~ norm(0,R)

%% Descrição do sistema discretizado

Ts = 1; % [s]

% Modelo verdadeiro
F = [ ...
    1 Ts 0; ...
    0 1 0; ...
    0 0 1 ...
    ];
H = [ 1 0 1 ];

q1 = 0.1^2; % Covariância do ruído de velocidade [m^2/s^2/Hz]
q2 = 0.1^2; % Covariância da incerteza do bias [m^2/Hz]

Qd = [ ...
    0 0 0; ...
    0 q1 0; ...
    0 0 q2 ...
    ]*Ts;
R = 0.025; % [m^2]

% Transformação linear modelo verdadeiro -> modelo de projeto
% a) Implementação ingênua

% Modelo de projeto (filtro)
F_linha = [ ...
    1 1; ...
    0 1 ...
    ];
H_linha = [ 1 0 ];
Qd_linha = [ ...
    q2*multiFactor 0; ...
    0 q1 ...
    ]*Ts;


%% Condições iniciais
P_posteriori_linha = [ ...
    1 0; ...
    0 0.25 ...
    ];

C_posteriori = [ ...
    0 0 0 0 0; ...
    0 0.25 0 0 0.25; ...
    0 0 1 1 0; ...
    0 0 1 1 0; ...
    0 0.25 0 0 0.25 ...
    ];


for i = 1:k
%% Equações da covariância computada pelo filtro
n = size(T,1);
% Obtemos P-_n a partir de P+_n-1
P_priori_linha = F_linha*P_posteriori_linha*F_linha'+Qd_linha;

% Obtemos K_n a partir de P-_n
K_linha = P_priori_linha*H_linha'/(H_linha*P_priori_linha*H_linha'+R);

% Obtemos P+_n a partir de P-_n
P_posteriori_linha = (eye(n)-K_linha*H_linha)*P_priori_linha;


%% Análise da covariância verdadeira 
n = size(T,2);

% Passo 1
DeltaF = T*F - F_linha*T;

% Passo 2
Fc = [F, zeros(length(DeltaF),length(F_linha)); DeltaF, F_linha];
Phic = Fc;

% Passo 3
Qc = [eye(n); T]*Qd*[eye(n), T'];

% Passo 4
DeltaH = H - H_linha*T;

% Passo 5
B = [...
    eye(n), zeros(n, length(K_linha*H_linha)); ...
    -K_linha*DeltaH, eye(length(K_linha*H_linha)) - K_linha*H_linha ...
    ];

% Passo 6
Kc = [ ...
    zeros(n,size(K_linha,2)); ...
    -K_linha ...
    ];

% Passo 7
C_priori = Phic*C_posteriori*Phic'+Qc;

% Passo 8
C_posteriori = B*C_priori*B' + Kc*R*Kc';

%% Coleta de resultados
% Propagada
resultados(1,i) = C_priori(4,4);
resultados(2,i) = C_priori(5,5);
resultados(3,i) = sqrt(P_posteriori_linha(1,1));
resultados(4,i) = sqrt(P_posteriori_linha(2,2));

% Atualizada
resultados(5,i) = C_posteriori(4,4);
resultados(6,i) = C_posteriori(5,5);
resultados(7,i) = sqrt(P_priori_linha(1,1));
resultados(8,i) = sqrt(P_priori_linha(2,2));
end

%% Plots
figure;
plot(resultados(1,:)', 'rx')
hold on
plot(resultados(3,:)', 'bx')
plot(resultados(5,:)', 'ro')
plot(resultados(7,:)', 'bo')
legend( ...
    'Segundo momento propagado', ...
    'Desvio padrão propagado', ...
    'Segundo momento atualizado', ...
    'Desvio padrão atualizado' ...
    );
title('Análise de covariância da posição');
grid on

figure;
plot(resultados(2,:)', 'rx')
hold on
plot(resultados(4,:)', 'bx')
plot(resultados(6,:)', 'ro')
plot(resultados(8,:)', 'bo')
legend( ...
    'Segundo momento propagado', ...
    'Desvio padrão propagado', ...
    'Segundo momento atualizado', ...
    'Desvio padrão atualizado' ...
    );
title('Análise de covariância da velocidade');
grid on

end

