%% Lista 1 - Parte 2
% Autor: Francisco Castro

ENTREGAVEL_parte1a;

%% Estabilização vertical
h_m = h;            % constante e igual à altitude inicial
estabVert.B = 0.01;
estabVert.C = 0.01;
estabVert.T_h = 60;

%% Modelamento da Terra e da gravidade
modTerra.R_0 = R_0;
modTerra.g_0 = g_0;
modTerra.e = e;

%% Integração numérica
y = [
    q0_base3;
    [0,0,0]';
    lambda;
    Lambda;
    h;
    h_m
    ];

%% Condições de parada
tf = 19.5*60;               % tempo final
i0 = 1;                     % index do início do movimento
i = i0;                     % index inicial
freqAquisicao = 100;        % Hz
p = 1/freqAquisicao;        % passo de integração
result = zeros(tf/p,length(y)+2);

%% Runge-kutta
for t=0:p:tf
    omega_B_i = w_b(i,:)';
    Asp_B_i = Asp_b(i,:)';
    
    k1 = f(t, y, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k2 = f(t + p/2, y + p/2*k1, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k3 = f(t + p/2, y + p/2*k2, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k4 = f(t + p, y + p*k3, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);

    y = y + p/6*(k1 + 2*k2 + 2*k3 + k4);
    result(i-i0+1,1) = t;       % tempo
    result(i-i0+1,2:12) = y;    % vetor de estados
    i = i+1;
end


%% Resultados
resultados.tempo = result(:,1)/60;
resultados.quaternion = result(:,2:5);
resultados.velocidade_NED = result(:,6:8);
resultados.latitude = result(:,9);
resultados.longitude = result(:,10);
resultados.altitude = result(:,11);

%% Plots
figure;
subplot(2,2,1)
plot(resultados.tempo, resultados.velocidade_NED)
legend('V_N','V_E','V_D');
title('Velocidade');
xlabel('Tempo [min]');

subplot(2,2,2)
plot(resultados.tempo, resultados.latitude);
title('Latitude');
xlabel('Tempo [min]');

subplot(2,2,3)
plot(resultados.tempo, resultados.longitude);
title('Longitude');
xlabel('Tempo [min]');

subplot(2,2,4)
plot(resultados.tempo, resultados.altitude);
title('Altitude');
xlabel('Tempo [min]');



%% TODO
% Estabilização canal vertical
% Vetor rotação ao final de cada volta com relação ao vetor original (q0)
% Índice de desalinhamento ao final de cada volta (escala log)
% Comparar índice de desalinhamento [arcseg] com magnitude do
% desalinhamento [arcseg] computado com quaternions
% Comparar os índices de normalidade e desalinhamento ao final de cada volta com e sem
% normalização de quaternion ao longo da integração (escala log)
% Converter as coordenadas geodésicas em cartesianas
% Apresentar a trajetória estimada pelo INS em relação ao ponto inicial

%3.1. Exiba o resultado da altura contra tempo, e o gráfico tri-dimensional da trajetória
%computada pelo INS.
% Plotar h x t
% Plotar 3d a trajetória

%3.2 Quantas voltas há no primeiro experimento? Plote Norte contra Leste na primeira volta
%e superponha sobre a imagem (do Google Maps) de satélite da montanha russa Montezum
%do parque de diversões HopiHari em Vinhedo, SP, nas duas primeiras voltas e até o final
%do experimento. O que ocorre?

%3.3. Qual o erro de posição relativa ao ponto inicial no plano horizontal e em altura ao
%final da primeira volta no primeiro experimento?

%3.4. Ainda com respeito ao primeiro experimento, o que ocorre entre o final da primeira
%volta e o início da segunda volta? E nas demais voltas, isto é, entre o fim de uma volta
%anterior e o início da seguinte?

%3.5. O que ocorre com a navegação inercial usando a IMU baseada em sensores MEMS do
%Xsens? 
% 

