resultados.tempo = result(:,1);
resultados.quaternion = result(:,2:5);
resultados.velocidade_NED = result(:,6:8);
resultados.latitude = result(:,9);
resultados.longitude = result(:,10);
resultados.altitude = result(:,11);

figure;

subplot(2,2,1)
plot(resultados.tempo/60, resultados.velocidade_NED)
legend('V_N','V_E','V_D');
title('Velocidade');
xlabel('Tempo [min]');

subplot(2,2,2)
plot(resultados.tempo/60, resultados.latitude);
title('Latitude');
xlabel('Tempo [min]');

subplot(2,2,3)
plot(resultados.tempo/60, resultados.longitude);
title('Longitude');
xlabel('Tempo [min]');

subplot(2,2,4)
plot(resultados.tempo/60, resultados.altitude);
title('Altitude');
xlabel('Tempo [min]');


%% Identificação das voltas
volta1.inicio = 0.8732;
volta1.fim = 3.201;
volta1.delta = 2.3278;

volta2.inicio = 4.212;
volta2.fim = 6.449;
volta2.delta = 2.237;

volta3.inicio = 7.552;
volta3.fim = 9.88;
volta3.delta = 2.328;

volta4.inicio = 10.92;
volta4.fim = 13.1;
volta4.delta = 2.18;
%%
index0 = 0;
correcaoPasso = 4;
indexFim.volta1 = round((19180 + index0)/correcaoPasso); % 4795
indexFim.volta2 = round((38680 + index0)/correcaoPasso); % 9670
indexFim.volta3 = round((59210 + index0)/correcaoPasso); % 14803
indexFim.volta4 = round((78910 + index0)/correcaoPasso); % 19728
%%
indexInicio.volta1 = 1263;
indexFim.volta1 = 4771;
indexFim.volta2 = 9615;
indexFim.volta3 = 14784;
indexFim.volta4 = 19654;
%%
indexFim.volta1 = 4663;
indexFim.volta2 = 9484;
indexFim.volta3 = 14646;
indexFim.volta4 = 19659;
%%
indexInicio.volta1 = 1261;
indexFim.volta1 = 4778;
indexFim.volta2 = 9632;
indexFim.volta3 = 14800;
indexFim.volta4 = 19800;

%%
indexFim
