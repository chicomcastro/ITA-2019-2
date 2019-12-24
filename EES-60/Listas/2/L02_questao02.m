clear

f = 1;
h = 1;
q = 0.01;
r = 1;

n = 1;
k = 50;

rng('default');
v = 0 + r.*randn(n,k);
w = 0 + q.*randn(n,k);
u = zeros(n,k)+1;
x(1) = 0;   % x0
z(1) = v(1);

result(1) = x(1);
for k = 2:50
    x(k) = f*x(k-1) + u(k-1) + w(k-1); % x_{k+1}
    z(k) = h*x(k) + v(k);
end

figure;
plot(1:k, x);
hold on;
plot(1:k, z);
grid on;
title("Compara��o");
xlabel("Tempo discretizado t");
ylabel("Estado x");
legend("Realiza��o", "Observa��o");

%%
clear

% Quantidade de itera��es
k = 50;
n = 1;

% Realiza��o
x = zeros(k,n);
u = (zeros(k,n)+1);
f = 1;
g = 1;

% Observa��o
z = zeros(k,n);
h = 1;

% Defini��o dos ru�dos
rng('default');
r = 1;
q = 0.01;
w = 0 + q.*randn(k,n);
v = 0 + r.*randn(k,n);

% Defini��o dos estimadores de m�nima covari�ncia
p00 = r/h;
P_posteriori = 10;      % P(0|0)

% Condi��es iniciais
x(1) = 0;
z(1) = v(1);

% Inicializa��o dos vetores de estima��o
x_priori = x;
x_posteriori = x;

% Evolu��o do estado
for i = 2:k
    % Realiza��o: o estado evolui para x_n
    x(i,:) = (f*x(i-1,:)' + g*u(i-1,:)' + w(i-1,:)')';
    
    % Estima��o: estimamos x-_n com base em x+_n-1
    x_priori(i,:) = (f*x_posteriori(i-1,:)'+g*u(i-1,:)')';
    
    % Obtemos P-_n a partir de P+_n-1
    P_priori = f*P_posteriori(i-1,:)*f'+q;
    
    % Obtemos K_n a partir de P-_n
    K = P_priori*h'/(h*P_priori*h'+r);
    
    % Obtemos P+_n a partir de P-_n
    P_posteriori(i,:) = (eye(n)-K*h)*P_priori;

    % Observa��o: medimos z_n
    z(i,:) = (h*x(i,:)' + v(i,:)')';
    
    % Previs�o (Filtro de Kalman): estimamos x+_n com base em x-_n, z_n e K_n 
    x_posteriori(i,:) = x_priori(i,:)+(K*(z(i,:)'-h*x_priori(i,:)'))';
    
end
%%
resultados(1:k,1) = 1:k;
resultados(1:k,end+1) = x;
resultados(1:k,end+1) = w;
resultados(1:k,end+1) = v;
resultados(1:k,end+1) = z;
resultados(1:k,end+1) = x_priori;
resultados(1:k,end+1) = P_priori;
resultados(1:k,end+1) = h*x_priori;
resultados(1:k,end+1) = h*P_priori*h'+r;
resultados(1:k,end+1) = z-h*x_priori;
resultados(1:k,end+1) = (z-h*x_priori)./((h*P_priori*h+r).^0.5);
resultados(1:k,end+1) = P_priori*h'*inv(h*P_priori*h'+r);
resultados(1:k,end+1) = x_posteriori;
resultados(1:k,end+1) = P_posteriori;
resultados(1:k,end+1) = x_posteriori./(P_posteriori.^0.5);
%%
figure;
plot(1:k, x);
hold on;
plot(1:k, x_posteriori);
hold on;
plot(1:k, z);
grid on;
xlabel("Tempo discretizado t");
ylabel("Estado x");
legend("Realiza��o", "Previs�o", "Observa��o");

figure
plot(1:k, x-x_posteriori, 'b');
hold on;
grid on;
plot(1:k, 3*sqrt(P_posteriori))
plot(1:k, -3*sqrt(P_posteriori))
xlabel("Tempo discretizado t");
ylabel("Erro absoluto de previs�o x(k) - x(k|k)");

