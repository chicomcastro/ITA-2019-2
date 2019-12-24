%Algoritmo Salichev
%Condições iniciais: V_NED(0), q_NED_b(0), ?(0), ?(0) ,h(0)

%% Preparação
tf = 19.5*60;               % tempo final

freqAquisicao = 100;        % Hz
T = 1/freqAquisicao;        % passo de integração
result = zeros(tf/(4*T),12);

% Pegando medida incrementais
w_b = w_b.*T;
Asp_b = Asp_b.*T;

%% Condições iniciais
q_NED_b = q0_base3;
V_NED = [0,0,0]';
lat = lambda;
lon = Lambda;
alt = h;
h_aux = h;

y = [...
    q_NED_b; ...
    V_NED; ...
    lat; ...
    lon; ...
    alt; ...
    h_aux ...
];

%% Integração
for k = 1:tf/(4*T)

    % Conjunto de 4 medidas
    omega_B_k = w_b(4*k-3 : 4*k-3+3, :)';
    Asp_B_k = Asp_b(4*k-3 : 4*k-3+3, :)';

    % Atualização de estado
    y = fs(k, y, omega_B_k, Asp_B_k, modTerra, estabVert, h_m, T);
    
    % Coleção de resultados
    result(k,1) = k;       % tempo
    result(k,2:12) = y;    % vetor de estados
    
end
