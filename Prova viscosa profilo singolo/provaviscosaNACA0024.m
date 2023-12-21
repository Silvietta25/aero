%% Grafico Cl(alpha) NACA 0024 simulazione viscosa
clear
close all
clc

% Si riportano i dati da XFoil
alpha=[0:0.5:11,12:0.5:25];
Cl= [0 0.0498 0.0990 0.1486 0.1978 0.2467 0.2956 0.3432 0.3913 0.4378 0.4842...
    0.5283 0.5695 0.6095 0.6516 0.6906 0.7294 0.7669 0.8113 0.8677 0.9330 0.9955...
    1.0567 1.1132 1.1108 1.1182 1.1313 1.1466 1.1615 1.1751 1.1869 1.1969 1.2050...
    1.2110 1.2152 1.2179 1.2196 1.2210 1.2219 1.2223 1.2234 1.2265 1.2284 1.2300...
    1.2343 1.2440 1.2441 1.2480 1.2598 1.2615];
plot(alpha, Cl, 'LineWidth',5,'Color','b')
hold on
plot(alpha,Cl(11)/alpha(11).*alpha,'LineWidth',5,'LineStyle','--','Color','r')
grid on
legend('Simulazione viscosa','Modello lineare','FontSize',40)
xlabel('alpha','FontSize',40)
ylabel('C_L','FontSize',40)
ax=gca;
ax.FontSize=40;
title('C_L(\alpha) NACA 0024 con simulazione viscosa XFoil','FontSize',50)