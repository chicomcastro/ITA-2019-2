%% Lista 01 - Parte 1b
% Autor: Francisco Castro
%%
% Nesta parte fazem-se as an�lises referentes aos dados presentes nos
% arquivos _SN500574Inside.mat_ e _SN500574Outside.mat_

%% Prepara��o do ambiente
close all
format
%clear
clc

%% Implementa��o da estrutura para o m�todo TRIAD
% Implementou-se, analogamente, o m�todo TRIAD considerando agora o vetor
% campo geomagn�tico ao inv�s do vetor velocidade angular para a correla��o
% entre as medidas no referencial do corpo e do NED. A metodologia �
% an�loga a feita para a parte 1a, por isso preferiu-se por omitir o
% c�digo.
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

%% Outside: Verifica��o de resultados e coment�rios
% Temos que os resultados expostos correspondem �s an�lises an�logas de
%
% # verifica��o dos determinantes das DCM iniciais obtidas
% # verifica��o das normas dos quaternions iniciais obtidos
% # desalinhamento para as matrizes ortogonais mais pr�ximas advindas do
% m�todo de ortonormaliza��o iterativa em rela��o a DCM original para cada base
% (desalinhamentoOrtoIt)
% # desalinhamento para as DCM provenientes dos quaternions obtidos como
% melhores quaternion de rota��o para cada estimativa de quaternion inicial
% (desalinhamentoNormIt) com rela��o � estimativa de DCM original para cada
% base
% # medida do desalinhamento cruzado entre as DCM originais, bem como as 
% obtidas do processo de
% ortonormaliza��o e normaliza��o iterativo comparadas com a estimativa de 
% DCM
% originalmente obtida utilizando a base 3, tomada como refer�ncia (melhor
% estimativa)
%analise;

%% Inside: Resultados
coord = coordOutside;
dados = data;
gm_NED = [17052.8,-6680.1,-13897.2];   % em nanoTesla (a ser normalizado)
non_incremental_metodologia;
%resultados;


%% Inside: Verifica��o de resultados e coment�rios
% Temos que os resultados expostos correspondem novamente 
% �s an�lises an�logas de
%
% # verifica��o dos determinantes das DCM iniciais obtidas
% # verifica��o das normas dos quaternions iniciais obtidos
% # desalinhamento para as matrizes ortogonais mais pr�ximas advindas do
% m�todo de ortonormaliza��o iterativa em rela��o a DCM original para cada base
% (desalinhamentoOrtoIt)
% # desalinhamento para as DCM provenientes dos quaternions obtidos como
% melhores quaternion de rota��o para cada estimativa de quaternion inicial
% (desalinhamentoNormIt) com rela��o � estimativa de DCM original para cada
% base
% # medida do desalinhamento cruzado entre as DCM originais, bem como as 
% obtidas do processo de
% ortonormaliza��o e normaliza��o iterativo comparadas com a estimativa de 
% DCM
% originalmente obtida utilizando a base 3, tomada como refer�ncia (melhor
% estimativa)

%%analise;

%%
% V�-se que os dados provenientes de dentro do laborat�rio parecem estar
% menos confi�veis (maiores desalinhamentos) quanto aos padr�es obtidos nas
% outras duas an�lises.

%% Conclus�es
% Assim como na Parte 1a, v�-se que as estimativas de DCM iniciais apresentam n�o 
% ortonormalidade associada � base que as gera, onde pode-se constatar pela
% compara��o entre as suas inversas e transpostas. Logo, uma base n�o
% ortogonal gera uma DCM tamb�m n�o ortonormal (esperado), o mesmo para
% uma base ortonormal.
%%
% Outra tend�ncia observada � que o desalinhamento (discrep�ncia) com
% rela��o as estimativas iniciais fica menor a partir da ortonormaliza��o
% iterativa a medida que se aumenta a complexidade da base utilizada para
% gerar a DCM. Tal medida ainda n�o nos diz qual a melhor escolha de base,
% por�m, considerando que a base 3 seria a ideal, uma vez que � ortonormal
% e conduz a uma DCM que � de fato uma matriz de rota��o e um quaternion
% que � de rota��o, temos que a discrep�ncia da base 2 para a base 3 �
% semelhante se comparada a da base 1, indicando ser a base 1 j� uma �tima
% candidata para uso em c�lculo sem a necessidade de apelar para mais
% elevados custos computacionais associados � utiliza��o das outras bases
% usando-se o m�todo de ortonormaliza��o iterativa.




