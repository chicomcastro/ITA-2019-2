clear

% Quantidades e dimens�es
k = 100;
n = 2;
T = 1;

% Realiza��o
x = zeros(k,n);
u = (zeros(k,n)+1);
F = [...
    1 T; ...
    0 1 ...
    ];
G = zeros(n);
TAU = [...
    T^2/2; ...
    T ...
    ];

% Observa��o
z = zeros(k,1);
H = [1 0];

% Defini��o dos ru�dos
rng('default');
mu = [ 0 0 ];
q = 0.01;
r = 1;
Q = q * eye(n);
R = r;
RQ = chol(Q);
RR = chol(R);
w = repmat(mu, k, n/2) + randn(k, n)*RQ;
v = randn(k, 1)*RR;

% Defini��o dos estimadores de m�nima covari�ncia
p00 = [...
    r r/T; ...
    r/T 2*r/T^2 ...
    ];
P_posteriori = p00;      % P(0|0)

% Condi��es iniciais
x(1,:) = [0 10];
z(1,:) = (H*x(1,:)' + v(1,:)')';

% Inicializa��o dos vetores de estima��o
x_priori = x;
x_posteriori = x;

% Evolu��o do estado
for i = 2:k
    % Realiza��o: o estado evolui para x_n
    x(i,:) = (F*x(i-1,:)' + G*u(i-1,:)' + w(i-1,:)')';
    
    % Estima��o: estimamos x-_n com base em x+_n-1
    x_priori(i,:) = (F*x_posteriori(i-1,:)'+G*u(i-1,:)')';
    
    % Obtemos P-_n a partir de P+_n-1
    P_priori = F*P_posteriori*F'+q;
    
    % Obtemos K_n a partir de P-_n
    K = P_priori*H'/(H*P_priori*H'+R);
    
    % Obtemos P+_n a partir de P-_n
    P_posteriori = (eye(n)-K*H)*P_priori;
    variancia(i,:) = [P_posteriori(1,1),P_posteriori(1,1)];

    % Observa��o: medimos z_n
    z(i,:) = (H*x(i,:)' + v(i,:)')';
    
    % Previs�o (Filtro de Kalman): estimamos x+_n com base em x-_n, z_n e K_n 
    x_posteriori(i,:) = x_priori(i,:)+(K*(z(i,:)'-H*x_priori(i,:)'))';
    
end

figure;
plot(1:k, x(:,1));
hold on;
plot(1:k, x_posteriori(:,1));
hold on;
plot(1:k, z(:,1));
grid on;
xlabel("Tempo discretizado k");
ylabel("Posi��o [m]");
legend("Realiza��o", "Previs�o", "Observa��o");

figure;
plot(1:k, x(:,2));
hold on;
plot(1:k, x_posteriori(:,2));
grid on;
xlabel("Tempo discretizado k");
ylabel("Velocidade [m/s]");
legend("Realiza��o", "Previs�o");

erro_medicao = x-x_posteriori;
variancia(1,:) = [p00(1,1),p00(2,2)];
figure
plot(1:k, erro_medicao(:,1), 'b');
hold on;
grid on;
plot(1:k, 3*sqrt(variancia(:,1)), 'r')
plot(1:k, -3*sqrt(variancia(:,1)), 'r')
title("Erro de previs�o de posi��o");
xlabel("Tempo discretizado k");
ylabel("Erro absoluto de previs�o x(k) - x(k|k) [m]");

figure
plot(1:k, erro_medicao(:,2), 'b');
hold on;
grid on;
plot(1:k, 3*sqrt(variancia(:,2)),'r')
plot(1:k, -3*sqrt(variancia(:,2)),'r')
title("Erro de previs�o de velocidade");
xlabel("Tempo discretizado k");
ylabel("Erro absoluto de previs�o v(k) - v(k|k) [m/s]");


