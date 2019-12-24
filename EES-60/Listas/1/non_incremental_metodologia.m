%% Metodologia para tratamento de dados não incrementais

%% Tratamento dos dados
% Os dados são herdados do script que chama este durante a sua execução,
% portanto, o tratamento aqui utiliza uma matriz de dados genérica.
t = dados(:,1);
w_b = dados(:,5:7);              % omega_bi, em rad/s
Asp_b = dados(:,2:4);            % força específica em m/s2
campo_geomag = dados(:,8:10);    % vet. campo geomagnético quase-normalizado


%% Condições iniciais
% Pegando a média das medidas dos 6 primeiros minutos (id <= 36000)
gravidade_ref = mean(sqrt(Asp_b(1:36000,1).^2 + Asp_b(1:36000,2).^2 + Asp_b(1:36000,3).^2));
gm = mean(campo_geomag(1:36000,:));
gm = gm/norm(gm);
gm_x = gm(1);
gm_y = gm(2);
gm_z = gm(3);
a_x = mean(Asp_b(1:36000,1));
a_y = mean(Asp_b(1:36000,2));
a_z = mean(Asp_b(1:36000,3));

%% Modelo de gravidade
% altere o parâmetro g0 no modelo de gravidade visto em sala de forma que a
% média da magnitude do vetor força específica nos primeiros seis minutos 
% seja igual à magnitude predita pelo modelo de gravidade nas respectivas 
% coordenadas iniciais
%%
referencia = gravidade_ref;
modeloGravidade;

%% Calibração do modelo de gravidade
lambda = coord.latitude;                % rad
h = coord.altitude;                     % m
gravidade_ref = referencia;
g_0 = double(subs(g_0));                % m/s2

%% Preparação para cálculo da DCM para o vetor campo geomagnético de referência
gm_NED = gm_NED/norm(gm_NED);
geomag_NED_x = gm_NED(1);
geomag_NED_y = gm_NED(2);
geomag_NED_z = gm_NED(3);