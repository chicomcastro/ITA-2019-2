%% Questão 1
load('resultados_montezum_salychev.mat')
load('modTerra.mat')

wgs84 = wgs84Ellipsoid('meter');
lat = resultados.latitude;
lon = resultados.longitude;
h = resultados.altitude;
[x,y,z] = geodetic2enu(lat,lon,h,lat(1),lon(1),h(1),wgs84,'radians');

%%
inicio = 1;
fim = indexInicio.volta1;
deltaT = resultados.tempo(fim) - resultados.tempo(inicio);
erroHorizontal = sqrt(x(fim)^2+y(fim)^2)
CASO3_estimativaDeriva(deltaT,Omega,modTerra,erroHorizontal)*180/pi*3600
CASO3_estimativaDeriva_SIMPLE(deltaT,modTerra,erroHorizontal)*180/pi*3600

%%
fim = indexFim.volta1;
deltaT = resultados.tempo(fim) - resultados.tempo(inicio);
erroHorizontal = sqrt(x(fim)^2+y(fim)^2)
CASO3_estimativaDeriva(deltaT,Omega,modTerra,erroHorizontal)*180/pi*3600
CASO3_estimativaDeriva_SIMPLE(deltaT,modTerra,erroHorizontal)*180/pi*3600

%%
fim = indexFim.volta4;
deltaT = resultados.tempo(fim) - resultados.tempo(inicio);
erroHorizontal = sqrt(x(fim)^2+y(fim)^2)
CASO3_estimativaDeriva(deltaT,Omega,modTerra,erroHorizontal)*180/pi*3600
CASO3_estimativaDeriva_SIMPLE(deltaT,modTerra,erroHorizontal)*180/pi*3600
