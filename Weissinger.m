clear
close all
clc

%% Dati iniziali

% velocità asintotica
U_inf = 1;

% angolo di sideslip (inserire in gradi)
beta = 0; 
beta = deg2rad(beta);

% Rotazione della velocità 
U_inf = [cos(beta) sin(beta) 0]' .* U_inf;

% Inserire numero di corpi
N_corpi = 1;

corda_radice = [1];
diedro = [0]; % inserire in gradi
diedro = deg2rad(diedro);

freccia = [20]; % inserire in gradi
freccia = deg2rad(freccia);

rastremazione = [0.8]; 
lambda = [10]; % allungamento alare
apertura_alare = 2;

LE_posizione_X = [0];
LE_posizione_Y = [0];
LE_posizione_Z = [0];

angolo_di_rollio = [0];
angolo_di_incidenza_aero = [0]; % alpha
angolo_di_beccheggio = [0];

% Discretizzazione
discretizzazione_semiapertura_alare = [20]; % numero panelli direzione apertura su SEMIALA
discretizzazione_corda = [10]; 

%% Parametri per costruzione ala

semi_apertura_alare = apertura_alare ./ 2;

superficie_alare = 2 * (semi_apertura_alare .* corda_radice .* ( 1 + rastremazione ) ./ 2);

superficie_piana = superficie_alare .* cos(diedro(1));

corda_tip = corda_radice .* rastremazione;

% Corda media aerodinamica MAC
MAC = (2/3) .* corda_radice .* ( (1 + rastremazione + rastremazione.^2)./(1 + rastremazione));

%% Creazione struttura ala

punti_di_controllo = [];
punti_di_induzione = [];
normali_pannelli = [];
vortici_a_infinito = [];
vortici = [];
mesh = [];
estremi_ala = [];

for i = 1:N_corpi
    Y_incidenza_aero = [cos(angolo_di_incidenza_aero(i))  0   sin(angolo_di_incidenza_aero(i));
                              0                      1               0                   ;          
                   -sin(angolo_di_incidenza_aero(i)) 0   cos(angolo_di_incidenza_aero(i))];
    X_diedro = [1 0 0;
                0 cos(diedro(i)) -sin(diedro(i));
                0 sin(diedro(i)) cos(diedro(i))];
    estremi_ala_radice_LE = [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]';
    estremi_ala = [estremi_ala, estremi_ala_radice_LE];

    estremi_ala_radice_TE = estremi_ala_radice_LE + Y_incidenza_aero * [corda_radice(i), 0, 0]';
    estremi_ala = [estremi_ala, estremi_ala_radice_TE];

    estremi_ala_tip_LE = [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]' + [0, semi_apertura_alare, 0]';
    estremi_ala_tip_LE = X_diedro * estremi_ala_tip_LE; % aggiunta diedro
    estremi_ala_tip_LE(1) = estremi_ala_radice_LE(1) + semi_apertura_alare .* tan(freccia(i)) + corda_radice - corda_tip; % aggiunta freccia
    estremi_ala_tip_LE = Y_incidenza_aero * estremi_ala_tip_LE;
    estremi_ala = [estremi_ala, estremi_ala_tip_LE];

    estremi_ala_tip_TE = estremi_ala_tip_LE + Y_incidenza_aero * [corda_tip 0 0]';

    estremi_ala = [estremi_ala, estremi_ala_tip_TE];

    % Discretizzazione direzione ala
    discretizzazione_direzione_ala = linspace(0, 1, discretizzazione_semiapertura_alare(i)+1);

    punti_LE = zeros(3, length(discretizzazione_direzione_ala));
    punti_TE = zeros(3, length(discretizzazione_direzione_ala));

    for j = 1:length(discretizzazione_direzione_ala)
        punti_LE(:, j) = (estremi_ala_tip_LE - estremi_ala_radice_LE).* discretizzazione_direzione_ala(j) + estremi_ala_radice_LE;
        punti_TE(:, j) = (estremi_ala_tip_TE - estremi_ala_radice_TE) .* discretizzazione_direzione_ala(j) + estremi_ala_radice_TE;
    end

    % Discretizzazione direzione corda
    discretizzazione_direzione_corda = linspace(0, 1, discretizzazione_corda(i)+1);

        punti_radice = zeros(3, length(discretizzazione_direzione_corda));
        punti_tip = zeros(3, length(discretizzazione_direzione_corda));

        quarto_tip = zeros(3, length(discretizzazione_direzione_ala));
        quarto_radice = zeros(3, length(discretizzazione_direzione_ala));

        trequarti_tip = zeros(3, length(discretizzazione_direzione_ala));
        trequarti_radice = zeros(3, length(discretizzazione_direzione_ala));

    for j = 1:length(discretizzazione_direzione_corda)
        punti_radice(:, j) = (estremi_ala_radice_TE - estremi_ala_radice_LE).* discretizzazione_direzione_corda(j) + estremi_ala_radice_LE;
        punti_tip(:, j) = (estremi_ala_tip_TE - estremi_ala_tip_LE) .* discretizzazione_direzione_corda(j) + estremi_ala_tip_LE;
    end

    for j = 1:length(discretizzazione_direzione_corda)-1
        quarto_tip(:, j) = (punti_tip(:, j) + 3.*punti_tip(:, j+1)) ./ 4;
        quarto_radice(:, j) = (punti_radice(:, j) + 3.*punti_radice(:, j+1)) ./ 4;
            
        trequarti_tip(:, j) = 3*(punti_tip(:, j) + 3.*punti_tip(:, j+1)) ./ 4;
        trequarti_radice(:, j) = 3*(punti_tip(:, j) + 3.*punti_tip(:, j+1)) ./ 4;
    end

    % Creazione mesh
    mesh = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 + 1);
    vortici = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 +1);
    punti_di_controllo = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    punti_di_induzione = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    normali_pannelli = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    vortici_a_infinito = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 + 1);
    
    lunghezza_vortici_infinito = 100 * corda_radice;
    
    % Discretizzazione direzione corda
    for j = 1:length(discretizzazione_direzione_corda)
    
        if j ~= length(discretizzazione_direzione_corda)
            for k = 1:length(discretizzazione_direzione_ala)
                vortici{j,k + length(discretizzazione_direzione_ala) - 1} = (quarto_tip(:, j) - quarto_radice(:, j)).*(discretizzazione_direzione_ala(k)) + quarto_radice(:, j);
                vortici_a_infinito{j,k + length(discretizzazione_direzione_ala) - 1} = vortici{j,k + length(discretizzazione_direzione_ala) - 1} + [lunghezza_vortici_infinito 0 0]';
    
                if k ~= length(discretizzazione_direzione_ala) 
                    punti_di_controllo{j,k + length(discretizzazione_direzione_ala) - 1} = (trequarti_tip(:, j) - trequarti_radice(:, j)).*(discretizzazione_direzione_ala(k)) + trequarti_radice(:, j);
                    punti_di_induzione{j,k + length(discretizzazione_direzione_ala) - 1} = (quarto_tip(:, j) - quarto_radice(:, j)).*((discretizzazione_direzione_ala(k) + discretizzazione_direzione_ala(k+1))./2) + quarto_radice(:, j);
                end
            end
        end

        for k = 1:length(discretizzazione_direzione_ala)
            mesh{j,k} = (punti_TE(k) - punti_LE(k)).* discretizzazione_direzione_corda(j) + punti_LE(k);
        end
    end

    % Creazione semiala per simmetria
    spost_X = LE_posizione_X(i);
    spost_Y = LE_posizione_Y(i);
    spost_Z = LE_posizione_Z(i);

    % Discretizzazione direzione corda
    for j = 1:length(discretizzazione_direzione_corda)

        if j ~= length(discretizzazione_direzione_corda)
            for k = 1:length(discretizzazione_direzione_ala)
                vortici{j,k} = -(vortici{j,2*length(discretizzazione_direzione_ala)-k} - [spost_X spost_Y spost_Z]') + [spost_X spost_Y spost_Z]';
                vortici_a_infinito{j,k} = -(vortici_a_infinito{j,2*length(discretizzazione_direzione_ala)-k} - [spost_X spost_Y spost_Z]') + [spost_X spost_Y spost_Z]';
    
                if k ~= length(discretizzazione_direzione_ala) 
                    punti_di_controllo{j,k} = -(punti_di_controllo{j,2*length(discretizzazione_direzione_ala)-1-k} - [spost_X spost_Y spost_Z]') + [spost_X spost_Y spost_Z]';
                    punti_di_induzione{j,k} = -(punti_di_induzione{j,2*length(discretizzazione_direzione_ala)-1-k} - [spost_X spost_Y spost_Z]') + [spost_X spost_Y spost_Z]';
                end
            end
        end
    end
end

figure;
grid on;
plot3(estremi_ala(1, :), estremi_ala(2, :), estremi_ala(3, :), '.');





