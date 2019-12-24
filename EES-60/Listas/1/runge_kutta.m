%% O método Runge–Kutta clássico de quarta ordem
%Um membro da família de métodos Runge–Kutta é usado com tanta frequência que costuma receber o nome de "RK4" ou simplesmente "o método Runge–Kutta".

% Seja um problema de valor inicial (PVI) especificado como segue:

%{\displaystyle y'=f(t,y),\quad y(t_{0})=y_{0}.} {\displaystyle y'=f(t,y),\quad y(t_{0})=y_{0}.}
%Então o método RK4 para este problema é dado pelas seguintes equações:

%{\displaystyle {\begin{aligned}y_{n+1}&=y_{n}+{h \over 6}\left(k_{1}+2k_{2}+2k_{3}+k_{4}\right)\\t_{n+1}&=t_{n}+h\\\end{aligned}}} {\displaystyle {\begin{aligned}y_{n+1}&=y_{n}+{h \over 6}\left(k_{1}+2k_{2}+2k_{3}+k_{4}\right)\\t_{n+1}&=t_{n}+h\\\end{aligned}}}
%onde {\displaystyle y_{n+1}} {\displaystyle y_{n+1}} é a aproximação por RK4 de {\displaystyle y(t_{n+1}),} {\displaystyle y(t_{n+1}),} e

%{\displaystyle {\begin{aligned}k_{1}&=f\left(t_{n},y_{n}\right)\\k_{2}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{1}\right)\\k_{3}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{2}\right)\\k_{4}&=f\left(t_{n}+h,y_{n}+hk_{3}\right)\\\end{aligned}}} {\displaystyle {\begin{aligned}k_{1}&=f\left(t_{n},y_{n}\right)\\k_{2}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{1}\right)\\k_{3}&=f\left(t_{n}+{h \over 2},y_{n}+{h \over 2}k_{2}\right)\\k_{4}&=f\left(t_{n}+h,y_{n}+hk_{3}\right)\\\end{aligned}}}
%Então, o próximo valor (yn+1) é determinado pelo valor atual (yn) somado com o produto do tamanho do intervalo (h) e uma inclinação estimada. A inclinação é uma média ponderada de inclinações:

%k1 é a inclinação no início do intervalo;
%k2 é a inclinação no ponto médio do intervalo, usando a inclinação k1 para determinar o valor de y no ponto tn + h/2 através do método de Euler;
%k3 é novamente a inclinação no ponto médio do intervalo, mas agora usando a inclinação k2 para determinar o valor de y;
%k4 é a inclinação no final do intervalo, com seu valor y determinado usando k3.

%%
%Ao fazer a média das quatro inclinações, um peso maior é dado para as inclinações no ponto médio:

%{\displaystyle {\mbox{inclinação}}={\frac {k_{1}+2k_{2}+2k_{3}+k_{4}}{6}}.} {\displaystyle {\mbox{inclinação}}={\frac {k_{1}+2k_{2}+2k_{3}+k_{4}}{6}}.}
%O método RK4 é um método de quarta ordem, significando que o erro por passo é da ordem de h5, enquanto o erro total acumulado tem ordem h4.

%% Implementação

%% Definições
f = @(t,y) sin(t);   % y' = f(t,y)
y = 1;          % y(t0) = y0
%t = 0;          % t0 = 0
h = 1e-3;       % passo

%% Condições de parada
n = 10e5;       % nro. iter
tf = 10;        % tempo final

i = 1;
result = zeros(tf/h,2);
%% Integração numérica
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