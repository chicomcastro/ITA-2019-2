%% O m�todo Runge�Kutta cl�ssico de quarta ordem
%Um membro da fam�lia de m�todos Runge�Kutta � usado com tanta frequ�ncia que costuma receber o nome de "RK4" ou simplesmente "o m�todo Runge�Kutta".

% Seja um problema de valor inicial (PVI) especificado como segue:

%{\displaystyle y'=f(t,y),\quad y(t_{0})=y_{0}.} {\displaystyle y'=f(t,y),\quad y(t_{0})=y_{0}.}
%Ent�o o m�todo RK4 para este problema � dado pelas seguintes equa��es:

%{\displaystyle {\begin{aligned}y_{n+1}&=y_{n}+{h \over 6}\left(k_{1}+2k_{2}+2k_{3}+k_{4}\right)\\t_{n+1}&=t_{n}+h\\\end{aligned}}} {\displaystyle {\begin{aligned}y_{n+1}&=y_{n}+{h \over 6}\left(k_{1}+2k_{2}+2k_{3}+k_{4}\right)\\t_{n+1}&=t_{n}+h\\\end{aligned}}}
%onde {\displaystyle y_{n+1}} {\displaystyle y_{n+1}} � a aproxima��o por RK4 de {\displaystyle y(t_{n+1}),} {\displaystyle y(t_{n+1}),} e

%{\displaystyle {\begin{aligned}k_{1}&=f\left(t_{n},y_{n}\right)\\k_{2}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{1}\right)\\k_{3}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{2}\right)\\k_{4}&=f\left(t_{n}+h,y_{n}+hk_{3}\right)\\\end{aligned}}} {\displaystyle {\begin{aligned}k_{1}&=f\left(t_{n},y_{n}\right)\\k_{2}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{1}\right)\\k_{3}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{2}\right)\\k_{4}&=f\left(t_{n}+h,y_{n}+hk_{3}\right)\\\end{aligned}}}
%Ent�o, o pr�ximo valor (yn+1) � determinado pelo valor atual (yn) somado com o produto do tamanho do intervalo (h) e uma inclina��o estimada. A inclina��o � uma m�dia ponderada de inclina��es:

%k1 � a inclina��o no in�cio do intervalo;
%k2 � a inclina��o no ponto m�dio do intervalo, usando a inclina��o k1 para determinar o valor de y no ponto tn + h/2 atrav�s do m�todo de Euler;
%k3 � novamente a inclina��o no ponto m�dio do intervalo, mas agora usando a inclina��o k2 para determinar o valor de y;
%k4 � a inclina��o no final do intervalo, com seu valor y determinado usando k3.

%%
%Ao fazer a m�dia das quatro inclina��es, um peso maior � dado para as inclina��es no ponto m�dio:

%{\displaystyle {\mbox{inclina��o}}={\frac {k_{1}+2k_{2}+2k_{3}+k_{4}}{6}}.} {\displaystyle {\mbox{inclina��o}}={\frac {k_{1}+2k_{2}+2k_{3}+k_{4}}{6}}.}
%O m�todo RK4 � um m�todo de quarta ordem, significando que o erro por passo � da ordem de h5, enquanto o erro total acumulado tem ordem h4.

%% Implementa��o

%% Defini��es
f = @(t,y) sin(t);   % y' = f(t,y)
y = 1;          % y(t0) = y0
%t = 0;          % t0 = 0
h = 1e-3;       % passo

%% Condi��es de parada
n = 10e5;       % nro. iter
tf = 10;        % tempo final

i = 1;
result = zeros(tf/h,2);
%% Integra��o num�rica
for t=0:h:tf
k1 = f(t,y);
k2 = f(t + h/2, y + h/2*k1);
k3 = f(t + h/2, y + h/2*k2);
k4 = f(t + h, y + h*k3);

y = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
result(i,1) = y;
result(i,2) = t;
i = i+1;
end

figure;
plot(result(:,2), result(:,1))