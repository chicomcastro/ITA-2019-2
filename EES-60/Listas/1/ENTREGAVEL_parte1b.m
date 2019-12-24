%% Lista 01 - Parte 1b
% Autor: Francisco Castro
%%
% Nesta parte fazem-se as análises referentes aos dados presentes nos
% arquivos _SN500574Inside.mat_ e _SN500574Outside.mat_

%% Preparação do ambiente
close all
format
%clear
clc

%% Implementação da estrutura para o método TRIAD
% Implementou-se, analogamente, o método TRIAD considerando agora o vetor
% campo geomagnético ao invés do vetor velocidade angular para a correlação
% entre as medidas no referencial do corpo e do NED. A metodologia é
% análoga a feita para a parte 1a, por isso preferiu-se por omitir o
% código.
triad_var;

%% Carregamento os dados
load('IMU Estacionaria/SN500574Inside.mat');
load('IMU Estacionaria/SN500574Outside.mat');

%% Coordenadas iniciais
coordInside.latitude = -(23.20995)*pi/180;
coordInside.longitude =  -(45.87712)*pi/180;
coordInside.altitude = 630; % m
coordOutside.latitude = -(23.20947)*pi/180;
coordOutside.longitude =  -(45.87722)*pi/180;
coordOutside.altitude = 630; % m

%% Outside: Resultados
coord = coordInside;
dados = outside;
gm_NED = [17046.4,-6679.6,-13904.3]';   % em nanoTesla (a ser normalizado)
non_incremental_metodologia;
resultados;

lambda = coordOutside.latitude;
Lambda = coordOutside.longitude;
h = coordOutside.altitude;

%% Outside: Verificação de resultados e comentários
% Temos que os resultados expostos correspondem às análises análogas de
%
% # verificação dos determinantes das DCM iniciais obtidas
% # verificação das normas dos quaternions iniciais obtidos
% # desalinhamento para as matrizes ortogonais mais próximas advindas do
% método de ortonormalização iterativa em relação a DCM original para cada base
% (desalinhamentoOrtoIt)
% # desalinhamento para as DCM provenientes dos quaternions obtidos como
% melhores quaternion de rotação para cada estimativa de quaternion inicial
% (desalinhamentoNormIt) com relação à estimativa de DCM original para cada
% base
% # medida do desalinhamento cruzado entre as DCM originais, bem como as 
% obtidas do processo de
% ortonormalização e normalização iterativo comparadas com a estimativa de 
% DCM
% originalmente obtida utilizando a base 3, tomada como referência (melhor
% estimativa)
%analise;

%% Inside: Resultados
coord = coordOutside;
dados = data;
gm_NED = [17052.8,-6680.1,-13897.2];   % em nanoTesla (a ser normalizado)
non_incremental_metodologia;
%resultados;


%% Inside: Verificação de resultados e comentários
% Temos que os resultados expostos correspondem novamente 
% às análises análogas de
%
% # verificação dos determinantes das DCM iniciais obtidas
% # verificação das normas dos quaternions iniciais obtidos
% # desalinhamento para as matrizes ortogonais mais próximas advindas do
% método de ortonormalização iterativa em relação a DCM original para cada base
% (desalinhamentoOrtoIt)
% # desalinhamento para as DCM provenientes dos quaternions obtidos como
% melhores quaternion de rotação para cada estimativa de quaternion inicial
% (desalinhamentoNormIt) com relação à estimativa de DCM original para cada
% base
% # medida do desalinhamento cruzado entre as DCM originais, bem como as 
% obtidas do processo de
% ortonormalização e normalização iterativo comparadas com a estimativa de 
% DCM
% originalmente obtida utilizando a base 3, tomada como referência (melhor
% estimativa)

%%analise;

%%
% Vê-se que os dados provenientes de dentro do laboratório parecem estar
% menos confiáveis (maiores desalinhamentos) quanto aos padrões obtidos nas
% outras duas análises.

%% Conclusões
% Assim como na Parte 1a, vê-se que as estimativas de DCM iniciais apresentam não 
% ortonormalidade associada à base que as gera, onde pode-se constatar pela
% comparação entre as suas inversas e transpostas. Logo, uma base não
% ortogonal gera uma DCM também não ortonormal (esperado), o mesmo para
% uma base ortonormal.
%%
% Outra tendência observada é que o desalinhamento (discrepância) com
% relação as estimativas iniciais fica menor a partir da ortonormalização
% iterativa a medida que se aumenta a complexidade da base utilizada para
% gerar a DCM. Tal medida ainda não nos diz qual a melhor escolha de base,
% porém, considerando que a base 3 seria a ideal, uma vez que é ortonormal
% e conduz a uma DCM que é de fato uma matriz de rotação e um quaternion
% que é de rotação, temos que a discrepância da base 2 para a base 3 é
% semelhante se comparada a da base 1, indicando ser a base 1 já uma ótima
% candidata para uso em cálculo sem a necessidade de apelar para mais
% elevados custos computacionais associados à utilização das outras bases
% usando-se o método de ortonormalização iterativa.




