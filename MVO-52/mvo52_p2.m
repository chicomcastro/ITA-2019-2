% MVO-52 - p2
clear
clc

dt = 1e-2;

l = 10;         % Comprimento do cabo do yoyo [m]
m = 1;          % Massa de 1U [kg]
s = 10e-2;      % Lado do cubesat [m]

Iu = m/6*s^2;   % Momento de inércia de 1U em torno de seu eixo [kg.m2]
A0 = 3*Iu;      % Momento de inércia de 3U em torno de x [kg.m2]
B0 = 3*Iu+2*m*s^2;% Momento de inércia de 3U em torno de y [kg.m2]
C0 = 3*Iu+2*m*s^2;% Momento de inércia de 3U em torno de z [kg.m2]

r0 = [ 0 s 0 ]';
w0 = [ 0 0 10 ]';
v0 = [ 0 0 0 ]';

A = A0;
B = B0;
C = C0;

I = [ ...
    A 0 0; ...
    0 B 0; ...
    0 0 C ...
    ];

r(1,:) = r0(:);
r(2,:) = -r0(:);
w = w0(:);
R = norm(r0);
vr(1,:) = v0;
vr(2,:) = -v0;
vt(1,:) = cross(w,r(1,:));
vt(2,:) = cross(w,r(2,:));

L = I*w + m*cross(r(1,:),vt(1,:))' + m*cross(r(2,:),vt(2,:))';
E = 1/2*C*norm(w)^2 + 1/2*m*norm(w)^2*norm(r(1,:))^2 + 1/2*m*norm(w)^2*norm(r(2,:))^2;

index = 1;

for i = 0:dt:10
    for j = 1:2
        %vt(j,:) = cross(w,r(j,:));
        r(j,:) = r(j,:) + vt(j,:)*dt;
        %r(j,:) = r(j,:)/norm(r(j,:))*R;

        if (norm(r(j,:)) < l)
            %vr(j,:) = vr(j,:) + norm(vt(j,:))^2/R*dt*r(j,:)/norm(r(j,:));
            %r(j,:) = r(j,:) + vr(j,:)*dt;
            %R = norm(r(j,:));
        end

        % Plota partícula
        particle(j).x = r(j,1);
        particle(j).y = r(j,2);
        particle(j).z = r(j,3);

        if (j == 1)
            color = 'b';
        else
            color = 'r';
        end
        
        plot(particle(j).x,particle(j).y, "o"+color);
        hold on;
    end
    
    % Atualiza momento de inércia do sistema
    A = 3*Iu + 2*m*(r(1,2)^2+r(1,3)^2);
    B = 3*Iu + 2*m*(r(1,1)^2+r(1,3)^2);
    C = 3*Iu + 2*m*(r(1,1)^2+r(1,2)^2);
    
    % Atualiza velocidade angular do sistema (só importa o eixo 3)
    w(3) = L(3)/(C+2*m*R^2);    %% DANDO PROBLEMA AQUI
    
    resultado(index,1:3) = vr(1,:);
    resultado(index,4:6) = w;
    index = index+1;
end

hold on;
plot(0,0, 'o');
axis([-2 2 -2 2]);
grid on;
axis equal;