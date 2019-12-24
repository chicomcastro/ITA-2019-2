%% Lista 01 - Parte 1a
% Autor: Francisco Castro
%%
% Nesta parte fazem-se as an�lises referentes aos dados presentes no
% arquivo _dados.dat_

%% Prepara��o do ambiente
close all
format
%clear
%clc

%% Implementa��o da estrutura para o m�todo TRIAD

%% 
% A express�o para a $DCM^{NED}_b$ utilizando o m�todo TRIAD pode ser
% encontrada considerando _inputs_ gen�ricos dos aceler�metros e gir�metros
% e montando os vetores em cada uma das bases no referencial
% do corpo e do NED que se correspondam para posterior associa��o e determina��o da respectiva DCM.
%%
% Portanto, temos
syms a_x a_y a_z w_x w_y w_z
Asp_b = [a_x ; a_y ; a_z];
Omega_b = [w_x; w_y; w_z];

syms g Omega lambda
Asp_NED = [0;0;-g];
Omega_NED = [Omega*cos(lambda); 0 ; -Omega*sin(lambda)];

%%
% Onde a base 1
%%
% $$[Asp, \Omega, Asp \times \Omega]$
%%
% que � uma base n�o ortonormal, pode ser definida da seguinte forma

A = [Asp_NED, Omega_NED, cross(Asp_NED,Omega_NED)];
B = [Asp_b, Omega_b, cross(Asp_b,Omega_b)];

D_NED_B_1 = B*inv(A);

%%
% Analogamente, temos a base 2
%%
% $$[Asp, Asp \times \Omega, Asp \times (Asp
% \times \Omega)]$
%%
% que � uma base ortogonal, mas n�o normalizada, descrita por
A = [Asp_NED, cross(Asp_NED,Omega_NED), cross(Asp_NED,cross(Asp_NED,Omega_NED))];
B = [Asp_b, cross(Asp_b,Omega_b), cross(Asp_b,cross(Asp_b,Omega_b))];

D_NED_B_2 = B*inv(A);

%%
% Ao passo que a base 3 fica 
%%
% $$\left [{ Asp \over |Asp|}, 
% {Asp \times \Omega \over |Asp \times \Omega|}, 
% {Asp \times (Asp \times \Omega) \over |Asp \times (Asp \times \Omega)|} \right ]$
%%
% que �, finalmente, uma base ortonormalizada.
A = [...
    Asp_NED/norm(Asp_NED), ...
    cross(Asp_NED,Omega_NED)/norm(cross(Asp_NED,Omega_NED)), ...
    cross(Asp_NED,cross(Asp_NED,Omega_NED))/norm(cross(Asp_NED,cross(Asp_NED,Omega_NED)))...
    ];
B = [...
    Asp_b/norm(Asp_b), ...
    cross(Asp_b,Omega_b)/norm(cross(Asp_b,Omega_b)), ...
    cross(Asp_b,cross(Asp_b,Omega_b))/norm(cross(Asp_b,cross(Asp_b,Omega_b)))...
    ];

D_NED_B_3 = B*inv(A);

%%
% Ao final dessa estrutura��o, temos as express�es simb�licas para a DCM
% NED -> B para as 3 bases propostas que poder�o ser calculadas ao final do
% modelamento.


%% Modelo de gravidade
% altere o par�metro g0 no modelo de gravidade visto em sala de forma que a
% m�dia da magnitude do vetor for�a espec�fica nos primeiros seis minutos 
% seja igual � magnitude predita pelo modelo de gravidade nas respectivas 
% coordenadas iniciais
%%
modeloGravidade;


%% Coordenadas iniciais
% Para as an�lises a serem feitas com os dados provenientes dos sensores
% instalados na montanha russa Montezum, ser�o consideradas as seguintes
% coordenadas
lambda = -(23 + 05/60 + 54.04/3600)*pi/180; % Latitude (rad)
Lambda =  -(47 + 00/60 + 41.55/3600)*pi/180; % Longitude (rad)
h = 774.6707;                               % Altitude (m)

%% Carregamento os dados
% Podemos, portanto, come�ar a mexer diretamente com os dados, aplicando-os
% nas express�es achadas para a DCM.

load('Montezum/dados.dat');

%%
% Tal que a frequ�ncia de aquisi��o de dados expressada pelas instru��es
% s�o tais que
frequenciaAmostragem = 100;    % Hz (frequ�ncia de amostragem)
dt = 1/frequenciaAmostragem;    % s  (intervalo de tempo entre medidas)

%% Tratamento dos dados
% Tendo os dados, podemos finalmente discrimin�-los quanto ao seu
% significado. Lembrando que os dados deste primeiro arquivo _dados.dat_
% apresentam medidas incrementais. Portanto, necessitou-se dividir pelo
% tempo entre amostragens _h_ para se obter a real grandeza a ser utilizada nas
% express�es obtidas.

id = dados(:,1);
tempo = id.*dt;                  % Tempo (s)
w_b = dados(:,2:4)./dt;      % Velocidade angular do corpo (rad/s)
Asp_b = dados(:,5:7)./dt;    % For�a espec�fica do corpo (m/s2)

%% Condi��es iniciais
% Pegando a m�dia das medidas dos 6 primeiros minutos (id <= 36000), temos
% os seguintes valores para a detemina��o da DCM na condi��o inicial.
%%
% Temos que � poss�vel fazer uma estimativa das grandezas de interesse em
% respeito ao NED a partir da norma da m�dia estacion�ria dos valores de for�a
% espec�fica e velocidade angular do corpo, dado que a sua correspond�ncia
% � expl�cita com a norma das mesmas grandezas no NED.
g = mean(sqrt(Asp_b(1:36000,1).^2 + Asp_b(1:36000,2).^2 + Asp_b(1:36000,3).^2));    %% M�dia da magnitude
Omega = norm(mean(w_b(1:36000,:)));     %% M�dia das medidas

%%
% Independentemente das estimativas acima, temos que as medidas a serem
% consideradas para velocidade angular e for�a espec�fica do corpo medidas
% pelos sensores no que diz respeito � determina��o da estimativa inicial
% da matriz DCM s�o, respectivamente

w_x = mean(w_b(1:36000,1));
w_y = mean(w_b(1:36000,2));
w_z = mean(w_b(1:36000,3));
a_x = mean(Asp_b(1:36000,1));
a_y = mean(Asp_b(1:36000,2));
a_z = mean(Asp_b(1:36000,3));

%%
% � importante ressaltar que, para um ru�do branco, a m�dia temporal dos 6
% primeiros minutos tende a reduzir a influ�ncia do mesmo, nos dando uma
% medida mais limpa para cada uma das grandezas que nos possibilita uma
% melhor determina��o da matriz de mudan�a de base que estamos
% interessados.


%% Calibra��o do modelo de gravidade
gravidade_ref = g;
g_0 = double(subs(g_0));

%% Resultados
% Dado todo o modelamento e tratamento realizados, temos, finalmente, como
% explicitar os valores atribu�dos para a velocidade angular e a for�a
% espec�fica do corpo nas matrizes DCM para cada uma das bases. Observe que
% os valores atribu�dos para as grandezas do corpo (assim como as do NED) s�o os mesmos para todas as bases, ou seja, a_x,
% a_y, a_z, w_x, w_y e w_z s�o invariantes (portanto, as estimativas para g e Omega tamb�m o s�o), o que muda � a formula��o dos
% vetores da base utilizada.

%%
% Com isso, temos que os resultados para a DCM inicial (D0), quaternion
% associado (q0) e �ngulos de Euler correspondentes, considerando
% rota��es 321 (os �ngulos est�o na ordem [yaw, pitch, roll] e s�o dados em radianos) s�o, para
% cada base.
%%
% *Base 1*
D0_NED_B_base1 = double(subs(D_NED_B_1));
q0_base1 = DCMtoQuaternion(D0_NED_B_base1);
euler0_base1 = quatToEuler(q0_base1);

%%
% *Base 2*
D0_NED_B_base2 = double(subs(D_NED_B_2));
q0_base2 = DCMtoQuaternion(D0_NED_B_base2);
euler0_base2 = quatToEuler(q0_base2);

%%
% *Base 3*
D0_NED_B_base3 = double(subs(D_NED_B_3));
q0_base3 = DCMtoQuaternion(D0_NED_B_base3);
euler0_base3 = quatToEuler(q0_base3);


%% Verifica��o de resultados
% A seguir s�o feitas algumas an�lises, verifica��es e avalia��es quanto �
% ortonormalidade das estimativas de D0 e a normalidade das estimativas dos
% respectivos quaternions.
%%
% *Ortonormalidade das DCMs*
%%
% Temos que a ortonormalidade das estimativas de DCM iniciais podem ser
% avaliadas atrav�s (i) do c�lculo do seu determinante (que deveria ser 1 para
% o caso de matrizes ortonormais) e (ii) a sua semelhan�a da sua matriz inversa
% com a a sua conjugada (defini��o de matriz ortonormal).
%%
% Primeiramente, temos que os determinantes das DCM obtidas s�o

determinantes = [det(D0_NED_B_base1), det(D0_NED_B_base2), det(D0_NED_B_base3)];

%%
% O que nos d� um primeiro ind�cio (i) de que o processo de ortonormaliza��o proposto pela
% redefini��o das bases consideradas realmente tem efeito e apresenta uma
% resposta, que pode ser mais ou menos desejada, ao custo computacional
% adicional empregado. Ademais, como � de se esperar, a matriz que n�o possui 
% determinante unit�rio n�o
% pode ser ortonormal e temos que, explicitamente, a que possui � (basta comparar inversa com transposta).

%%
% *Normalidade*
%%
% Quanto � normalidade dos quaternions, temos que as normas avaliadas para
% cada base s�o
normasQuat = [norm(q0_base1), norm(q0_base2), norm(q0_base3)];

%%
% Que nos confere um resultado an�logo ao anterior: apenas a �ltima � um
% quaternion de rota��o, logo, h� um pre�o a ser pago por se utilizar cada
% um dos dois primeiros como um quaternion de rota��o para se fazer os c�lculos.
% Tal pre�o pode ser contabilizado de algumas formas, como
% por exemplo, pela norma da diferen�a entre o quaternion e o quaternion de
% rota��o mais pr�ximo obtido por um processo de normaliza��o iterativa.
% Esta an�lise ser� afeita adiante.

%% An�lises de normaliza��o e ortonormaliza��o iterativa
% Em seguida mostra-se os resultados para o m�todo
% de ortonormaliza��o iterativo para obter a matriz ortonormal mais pr�xima
% das estimativas da DCM computadas pelo TRIAD, bem como os resultaados
% an�logos para o m�todo de normaliza��o iterativo para obter o quaternion 
% de rota��o mais pr�ximo da estimativa computada via TRIAD.

%%
% *Ortonormaliza��o iterativa*
%%
% Temos que as matrizes ortonormais mais pr�ximas �s obtidas apresentam,
% conforme esperado, as propriedades de uma matriz ortonormal. Mas
% analisar-se-� o seu desalinhamento com rela��o � DCM que a originou para
% tentar metrificar a magnitude do erro na orienta��o resultante do uso de
% uma em detrimento da outra. Dessa forma, temos
D0_base1_orto = ortonormalizacaoIterativo(D0_NED_B_base1);
D0_base2_orto = ortonormalizacaoIterativo(D0_NED_B_base2);
D0_base3_orto = ortonormalizacaoIterativo(D0_NED_B_base3);
format shortg;

desalinhamentoOrtoIt = [...
    desalinhamento(D0_base1_orto, D0_NED_B_base1),...
    desalinhamento(D0_base2_orto, D0_NED_B_base2),...
    desalinhamento(D0_base3_orto, D0_NED_B_base3) ...
    ];
%%
% Onde v�-se que a normaliza��o iterativa causa diferencia��o similar (da
% mesma ordem de grandeza) quanto a matriz original independentemente se a 
% base utilizada � ortogonal (base 2) ou n�o (base 2). Tal m�trica n�o nos
% permite invalidar qualquer uma das 3 possibilidades, mas nos d� ind�cios
% de que a ortonormaliza��o iterativa n�o nos d� tanto ganho a mais de
% precis�o mesmo utilizando uma base que necessita de maior poder
% computacional.
    
%%
% *Normaliza��o iterativa*
%%
% Quanto ao tratamento dos quaternions, optou-se por normalizar os
% resultados para estimativa de quaternion inicial em cada caso e comparar,
% adotando como m�trica o desalinhamento do novo quaternion de rota��o obtido
% com o quaternion original. Com isso, teremos
% uma compara��o an�loga � feita acima e poderemos avaliar, a grosso modo,
% os benef�cios de uma em rela��o �s outras quanto aos quaternions gerados.
%%
q_rot1 = normalizacaoIterativa(q0_base1);
q_rot2 = normalizacaoIterativa(q0_base2);
q_rot3 = normalizacaoIterativa(q0_base3);

desalinhamentoNormIt = [...
    desalinhamento(quatToDCM(q_rot1), D0_NED_B_base1),... % Calcula o desalinhamento
    desalinhamento(quatToDCM(q_rot2), D0_NED_B_base2),... % Calcula o desalinhamento
    desalinhamento(quatToDCM(q_rot3), D0_NED_B_base3) ... % Calcula o desalinhamento
    ];
%%
% Que nos d� um resultado an�logo ao anterior: h� pouco ganho, em termos de
% diferen�a do resultado original (medida da discrep�ncia) por se usar uma
% base computacionalmente mais custosa de se calcular. Curiosamente a
% varia��o de orienta��o (desalinhamento com rela��o ao uso da DCM proveniente do
% quaternion de rota��o mais pr�ximo ao inv�s da DCM original)
% que ter-se-� por normalizar iterativamente o quaternion � praticamente a
% mesma, considerando a m�trica adotada, para a base 1 e 2.
%%
% Outra an�lise que pode ser feita � o desalinhamento cruzado entre os
% resultados para cada base comparativa aos resultados originais da base 3,
% que pode ser considerada como a melhor estimativa a ser utilizada para os
% c�lculos necess�rios. Nesse sentido, temos que
desalinhamentoCruzado = [...
    % Linha 1: compara��o estimativas originais com a base 3
    desalinhamento(D0_NED_B_base1, D0_NED_B_base3),...
    desalinhamento(D0_NED_B_base2, D0_NED_B_base3),...
    desalinhamento(D0_NED_B_base3, D0_NED_B_base3);...
    % Linha 2: compara��o matriz ortonormal mais pr�xima com a base 3
    desalinhamento(D0_base1_orto, D0_NED_B_base3),...
    desalinhamento(D0_base2_orto, D0_NED_B_base3),...
    desalinhamento(D0_base3_orto, D0_NED_B_base3);...
    % Linha 3: compara��o DCM do quaternion de rota��o mais pr�ximo c/ a base 3
    desalinhamento(quatToDCM(q_rot1), D0_NED_B_base3),...
    desalinhamento(quatToDCM(q_rot2), D0_NED_B_base3),...
    desalinhamento(quatToDCM(q_rot3), D0_NED_B_base3)...    
    ];
%%
% Que nos fornece uma an�lise mais interessante: os processos de ortonormaliza��o
% iterativa das DCMs resultantes da base 1 e 2 nos d� praticamente o
% que obter�amos fazendo os c�lculos e considerando os resultados originais
% para a base 3. Ou seja, basta, nesse sentido, avaliar qual dessas
% alternativas possui menor custo computacional e, � grosso modo, a 
% precis�o das duas ser�
% semelhante se usada �lgebra matricial (DCMs). Caso opte-se por usar
% �lgebra de quaternions, que � menos custosa, h� ainda um erro associado
% que � consider�vel se usados os resultados da normaliza��o para as bases
% 1 e 2.
