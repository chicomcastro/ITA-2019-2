%% Lista 01 - Parte 1b
% Autor: Francisco Castro

%% Preparação do ambiente
clear
clc

%% Implementação da estrutura para o método TRIAD
triad_var

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

%% Outside
coord = coordInside;
dados = outside;
gm_NED = [17046.4,-6679.6,-13904.3]';
non_incremental_metodologia;
resultados;

%% Inside
coord = coordOutside;
dados = data;
gm_NED = [17052.8,-6680.1,-13897.2];
non_incremental_metodologia;
resultados;



