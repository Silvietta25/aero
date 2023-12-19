close all
clear all
clc

%% Input
U_Inf_Mag = 1;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;
config.NCorpi = 2; % numero corpi

%% Dati del tandem
config.RootChord = [1 0.5]; % corde delle 2 ali
config.DihedralAngle = [0 0]; % angolo di diedro in gradi delle 2 ali
config.SweepAngle = [10 10]; % angolo di freccia in gradi delle 2 ali
config.TaperRatio = [1 0.5]; % rapporto di rastremazione delle 2 ali
config.Span = [2 1]; % apertura alare delle 2 ali
config.LEPosition_X = [0 5]; % posizione su x dell'origine delle 2 ali
config.LEPosition_Y = [0 0]; % posizione su y dell'origine delle 2 ali
config.LEPosition_Z = [0 1]; % posizione su z dell'origine delle 2 ali

config.RotationAngle_X = [0 0]; % angolo di rollio (lo lascio sempre nullo)
config.RotationAngle_Z = [0 0]; % angolo d'imbardata (lo lascio sempre nullo)

% Discretization options
config.SemiSpanwiseDiscr = [5 5]; % discretizzazione su semiapertura alare per le 2 ali
config.ChordwiseDiscr = [5 5]; % discretizzazione su corda per le 2 ali

%% Salvo dati dei coefficienti aerodinamici
% Cl1(i,j) sarà il Cl del corpo 1 con incidenza ALPHA1=ALPHA1_vect(i) e con
% ALPHA2=ALPHA2_vect(j)
ALPHA1_vect=4:0.5:10;
ALPHA2_vect=20:0.5:25;
Cl1=zeros(length(ALPHA1_vect),length(ALPHA2_vect));
Cl2=zeros(length(ALPHA1_vect),length(ALPHA2_vect));
Cd1=zeros(length(ALPHA1_vect),length(ALPHA2_vect));
Cd2=zeros(length(ALPHA1_vect),length(ALPHA2_vect));
for i=1:length(ALPHA1_vect)
    ALPHA1=ALPHA1_vect(i);
    for j=1:length(ALPHA2_vect)
        ALPHA2=ALPHA2_vect(j);
        [Cl_3D_Corpo1,Cl_3D_Corpo2,Cd_3D_Corpo1,Cd_3D_Corpo2]=CL_e_CD_tandem...
        (config,ALPHA1,ALPHA2,U_Inf,rho);
        Cl1(i,j)=Cl_3D_Corpo1;
        Cl2(i,j)=Cl_3D_Corpo2;
        Cd1(i,j)=Cd_3D_Corpo1;
        Cd2(i,j)=Cd_3D_Corpo2;
    end
end

%% Grafico con alpha1 variabile e alpha2 fisso ad esempio
% Rappresento, in pratica, una colonna di Cl corrispondente all'alpha2
% voluto (esempio: prendo la prima colonna, cioè alpha2=2)--> ho il grafico
% di Cl e Cd in funzione di ALPHA1_vect, con ALPHA2 fissato
figure
plot(ALPHA1_vect,Cl1(:,1)','LineWidth',2)
hold on
plot(ALPHA1_vect,Cl2(:,1)','LineWidth',2)
hold on
plot(ALPHA1_vect,Cd1(:,1)','LineWidth',2)
hold on
plot(ALPHA1_vect,Cd2(:,1)','LineWidth',2)
legend('Cl1','Cl2','Cd1','Cd2')
grid on
title('Coefficienti aerodinamici Weissinger tandem con \alpha_2 fissato')

%% Grafico con alpha2 variabile e alpha1 fisso ad esempio
% Rappresento, in pratica, una riga di Cl corrispondente all'alpha1
% voluto (esempio: prendo la prima riga, cioè alpha1=4)--> ho il grafico
% di Cl e Cd in funzione di ALPHA2_vect, con ALPHA1 fissato
figure
plot(ALPHA2_vect,Cl1(1,:)','LineWidth',2)
hold on
plot(ALPHA2_vect,Cl2(1,:)','LineWidth',2)
hold on
plot(ALPHA2_vect,Cd1(1,:)','LineWidth',2)
hold on
plot(ALPHA2_vect,Cd2(1,:)','LineWidth',2)
legend('Cl1','Cl2','Cd1','Cd2')
grid on
title('Coefficienti aerodinamici Weissinger tandem con \alpha_1 fissato')


