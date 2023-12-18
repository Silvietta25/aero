%% CODICE WEISSINGER

clear
close all
clc

%% Dati iniziali

% velocità asintotica
U_inf = [1 0 0]';

% angolo di sideslip (inserire in gradi)
beta = 0; 
beta = deg2rad(beta);

% Rotazione della velocità 
U_inf = [cos(beta) sin(beta) 0]' .* U_inf;

% Inserire numero di corpi
N_corpi = 2;

corda_radice = [1 0.5]';
diedro = [0 0]'; % inserire in gradi
diedro = deg2rad(diedro);

freccia = [10 10]'; % inserire in gradi
freccia = deg2rad(freccia);

rastremazione = [1 0.5]'; 
apertura_alare = [2 1]';

LE_posizione_X = [0 5]';
LE_posizione_Y = [0 0]';
LE_posizione_Z = [0 1]';

% angolo_di_rollio = [0];
angolo_di_incidenza_aero = [10 25]'; % alpha
% angolo_di_beccheggio = [0];

% angolo_di_rollio = deg2rad(angolo_di_rollio);
angolo_di_incidenza_aero = deg2rad(angolo_di_incidenza_aero); % alpha
% angolo_di_beccheggio =deg2rad(angolo_di_beccheggio);

% Discretizzazione
discretizzazione_semiapertura_alare = [5 5]'; % numero panelli direzione apertura su SEMIALA
discretizzazione_corda = [5 5]'; 

% Inizializzazione sistema lineare
num_tot_pann = 0;
for i = 1:N_corpi
    num_tot_pann = num_tot_pann + 2*discretizzazione_corda(i)*discretizzazione_semiapertura_alare(i);
end
A = zeros(num_tot_pann, num_tot_pann);
b = zeros(num_tot_pann, 1);

gamma = cell(N_corpi,1);


%% Creazione struttura ala
estremi_ala = cell(1,N_corpi);

mesh = cell(1,N_corpi);
vortici = cell(1,N_corpi);
punti_di_controllo = cell(1,N_corpi);
punti_di_induzione = cell(1,N_corpi);
normali_pannelli = cell(1,N_corpi);
vortici_a_infinito = cell(1,N_corpi);

figure

for i = 1:N_corpi

    %% CREAZIONE GEOMETRIA ALA
    % Parametri per costruzione ala
    
    semi_apertura_alare = apertura_alare(i) / 2;
    
    superficie_alare = 2 * (semi_apertura_alare * corda_radice(i) * ( 1 + rastremazione(i) ) / 2);
    
    superficie_piana = superficie_alare * cos(diedro(i));
    
    corda_tip = corda_radice(i) * rastremazione(i);
    
    % Corda media aerodinamica MAC
    MAC = (2/3) * corda_radice(i) * ( (1 + rastremazione(i) + rastremazione(i)^2)/(1 + rastremazione(i)));

    Y_incidenza_aero = [cos(angolo_di_incidenza_aero(i))        0   sin(angolo_di_incidenza_aero(i));
                                        0                       1                   0               ;          
                        -sin(angolo_di_incidenza_aero(i))       0   cos(angolo_di_incidenza_aero(i))];

    X_diedro = [    1              0                  0      ;
                    0       cos(diedro(i))    -sin(diedro(i));
                    0       sin(diedro(i))    cos(diedro(i))];

    estremi_ala_radice_LE = [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]';
    estremi_ala{1,i} = [estremi_ala{1,i}, estremi_ala_radice_LE];

    estremi_ala_radice_TE = estremi_ala_radice_LE + Y_incidenza_aero * [corda_radice(i), 0, 0]';
    estremi_ala{1,i} = [estremi_ala{1,i}, estremi_ala_radice_TE];

    estremi_ala_tip_LE = [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]' + [0, semi_apertura_alare, 0]';
    estremi_ala_tip_LE = X_diedro * estremi_ala_tip_LE; % aggiunta diedro
    estremi_ala_tip_LE(1) = estremi_ala_radice_LE(1) + semi_apertura_alare * tan(freccia(i)) + corda_radice(i)/4 - corda_tip/4; % aggiunta freccia
    estremi_ala_tip_LE = Y_incidenza_aero * (estremi_ala_tip_LE - [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]') + [LE_posizione_X(i) LE_posizione_Y(i) LE_posizione_Z(i)]';
    estremi_ala{1,i} = [estremi_ala{1,i}, estremi_ala_tip_LE];

    estremi_ala_tip_TE = estremi_ala_tip_LE + Y_incidenza_aero * [corda_tip 0 0]';

    estremi_ala{1,i} = [estremi_ala{1,i}, estremi_ala_tip_TE];

    % Discretizzazione direzione ala
    discretizzazione_direzione_ala = linspace(0, 1, discretizzazione_semiapertura_alare(i)+1);

    punti_LE = zeros(3, length(discretizzazione_direzione_ala));
    punti_TE = zeros(3, length(discretizzazione_direzione_ala));

    for j = 1:length(discretizzazione_direzione_ala)
        punti_LE(:, j) = (estremi_ala_tip_LE - estremi_ala_radice_LE)* discretizzazione_direzione_ala(j) + estremi_ala_radice_LE;
        punti_TE(:, j) = (estremi_ala_tip_TE - estremi_ala_radice_TE) * discretizzazione_direzione_ala(j) + estremi_ala_radice_TE;
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
        quarto_tip(:, j) = (punti_tip(:, j+1) + 3.*punti_tip(:, j)) ./ 4;
        quarto_radice(:, j) = (punti_radice(:, j+1) + 3.*punti_radice(:, j)) ./ 4;
            
        trequarti_tip(:, j) = 3*((punti_tip(:, j+1) - punti_tip(:, j)) ./ 4) + punti_tip(:, j);
        trequarti_radice(:, j) = 3*((punti_radice(:, j+1) - punti_radice(:, j)) ./ 4) + punti_radice(:, j);
    end

    % Creazione mesh
    mesh_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 + 1);
    vortici_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 +1);
    punti_di_controllo_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    punti_di_induzione_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    normali_pannelli_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2);
    vortici_a_infinito_i = cell(discretizzazione_corda(i), discretizzazione_semiapertura_alare(i)*2 + 1);
    
    lunghezza_vortici_infinito = 100 * corda_radice;
    
    % Discretizzazione direzione corda
    for j = 1:length(discretizzazione_direzione_corda)
    
        if j ~= length(discretizzazione_direzione_corda)
            for k = 1:length(discretizzazione_direzione_ala)
                vortici_i{j,k + length(discretizzazione_direzione_ala) - 1} = (quarto_tip(:, j) - quarto_radice(:, j)).*(discretizzazione_direzione_ala(k)) + quarto_radice(:, j);
                vortici_a_infinito_i{j,k + length(discretizzazione_direzione_ala) - 1} = vortici_i{j,k + length(discretizzazione_direzione_ala) - 1} + [lunghezza_vortici_infinito(i) 0 0]';
    
                if k ~= length(discretizzazione_direzione_ala) 
                    punti_di_controllo_i{j,k + length(discretizzazione_direzione_ala) - 1} = (trequarti_tip(:, j) - trequarti_radice(:, j)).*((discretizzazione_direzione_ala(k)+ discretizzazione_direzione_ala(k+1))./2) + trequarti_radice(:, j);
                    punti_di_induzione_i{j,k + length(discretizzazione_direzione_ala) - 1} = (quarto_tip(:, j) - quarto_radice(:, j)).*((discretizzazione_direzione_ala(k) + discretizzazione_direzione_ala(k+1))./2) + quarto_radice(:, j);
                end
            end
        end

        for k = 1:length(discretizzazione_direzione_ala)
            mesh_i{j,k + length(discretizzazione_direzione_ala) - 1} = (punti_TE(:,k) - punti_LE(:,k)).* discretizzazione_direzione_corda(j) + punti_LE(:,k);
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
                vortici_i{j,k} = (vortici_i{j,2*length(discretizzazione_direzione_ala)-k} - [spost_X spost_Y spost_Z]');
                vortici_i{j,k}(2) = -(vortici_i{j,k}(2));
                vortici_i{j,k} = vortici_i{j,k} + [spost_X spost_Y spost_Z]';
                vortici_a_infinito_i{j,k} = (vortici_a_infinito_i{j,2*length(discretizzazione_direzione_ala)-k} - [spost_X spost_Y spost_Z]');
                vortici_a_infinito_i{j,k}(2) = -(vortici_a_infinito_i{j,k}(2));
                vortici_a_infinito_i{j,k} = vortici_a_infinito_i{j,k} + [spost_X spost_Y spost_Z]';
    
                if k ~= length(discretizzazione_direzione_ala) 
                    punti_di_controllo_i{j,k} = (punti_di_controllo_i{j,2*length(discretizzazione_direzione_ala)-1-k} - [spost_X spost_Y spost_Z]');
                    punti_di_controllo_i{j,k}(2) = -(punti_di_controllo_i{j,k}(2));
                    punti_di_controllo_i{j,k} = punti_di_controllo_i{j,k} + [spost_X spost_Y spost_Z]';
                    punti_di_induzione_i{j,k} = (punti_di_induzione_i{j,2*length(discretizzazione_direzione_ala)-1-k} - [spost_X spost_Y spost_Z]');
                    punti_di_induzione_i{j,k}(2) = -(punti_di_induzione_i{j,k}(2));
                    punti_di_induzione_i{j,k} = punti_di_induzione_i{j,k} + [spost_X spost_Y spost_Z]';
                end
            end
        end

        for k = 1:length(discretizzazione_direzione_ala) - 1
            mesh_i{j,k} = (mesh_i{j,2*length(discretizzazione_direzione_ala)-k} - [spost_X spost_Y spost_Z]');
            mesh_i{j,k}(2) = -(mesh_i{j,k}(2));
            mesh_i{j,k} = mesh_i{j,k} + [spost_X spost_Y spost_Z]';
        end
    end
    
    % plot geometria
    pannelli_x = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    pannelli_y = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    pannelli_z = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    
    vortici_x = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    vortici_y = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    vortici_z = zeros(length(discretizzazione_direzione_corda), 2*length(discretizzazione_direzione_ala)-1);
    
    for j = 1:length(discretizzazione_direzione_corda)
        for k = 1:2*length(discretizzazione_direzione_ala)-1
            pannelli_x(j, k) = mesh_i{j,k}(1);
            pannelli_y(j, k) = mesh_i{j,k}(2);
            pannelli_z(j, k) = mesh_i{j,k}(3);
            
            if j ~= length(discretizzazione_direzione_corda)
                vortici_x(j, k) = vortici_i{j,k}(1);
                vortici_y(j, k) = vortici_i{j,k}(2);
                vortici_z(j, k) = vortici_i{j,k}(3);
            end
        end
    end

    for j = 1:length(discretizzazione_direzione_corda)
        plot3(pannelli_x(j, :), pannelli_y(j, :), pannelli_z(j, :), 'color', 'k');
        hold on
        plot3(vortici_x(j, :), vortici_y(j, :), vortici_z(j, :), 'color', "#4DBEEE");
        hold on
    end
    for k = 1:2*length(discretizzazione_direzione_ala)-1
        plot3(pannelli_x(:, k), pannelli_y(:, k), pannelli_z(:, k), 'Color', 'k');
        hold on
    end
    grid on
    axis equal

    for j = 1:discretizzazione_corda(i)
        for k = 1:2*discretizzazione_semiapertura_alare(i)
            plot3(punti_di_controllo_i{j,k}(1), punti_di_controllo_i{j,k}(2), punti_di_controllo_i{j,k}(3), '.' ,'color', 'red')
            plot3(punti_di_induzione_i{j,k}(1), punti_di_induzione_i{j,k}(2), punti_di_induzione_i{j,k}(3), '.' ,'color', 'blue')
        end
    end

    % Calcolo normali pannelli
    for j = 1:discretizzazione_corda(i)
        for k = 1:2*discretizzazione_semiapertura_alare(i)
            normali_pannelli_i{j,k} = cross(mesh_i{j,k}-mesh_i{j+1,k+1}, mesh_i{j,k+1}-mesh_i{j+1,k});
            normali_pannelli_i{j,k} = -normali_pannelli_i{j,k}/norm(normali_pannelli_i{j,k});
        end
    end

    mesh{1,i} = mesh_i;
    vortici{1,i} = vortici_i;
    punti_di_controllo{1,i} = punti_di_controllo_i;
    punti_di_induzione{1,i} = punti_di_induzione_i;
    normali_pannelli{1,i} = normali_pannelli_i;
    vortici_a_infinito{1,i} = vortici_a_infinito_i;
    
end


%% SISTEMA LINEARE
% Costruzione matrice
pannello_indotto = 0;

for i = 1:N_corpi
    for j_indotto = 1:discretizzazione_corda(i)
        for k_indotto = 1:2*discretizzazione_semiapertura_alare(i)

            pannello_indotto = pannello_indotto + 1;
            vortice_inducente = 0;

            punto_controllo = punti_di_controllo{1,i}{j_indotto, k_indotto};
            normale= normali_pannelli{1,i}{j_indotto, k_indotto};

            
            for i_inducente = 1:N_corpi
                for j_inducente = 1:discretizzazione_corda(i_inducente)
                    for k_inducente = 1:2*discretizzazione_semiapertura_alare(i_inducente)

                        vortice_inducente = vortice_inducente + 1;
                        U_ind = 0;

                        % Induzione vortice semi-infinito sx
                        estremo_in = vortici{1,i_inducente}{j_inducente, k_inducente};
                        estremo_fin = vortici_a_infinito{1,i_inducente}{j_inducente, k_inducente};
                        r0 = estremo_in - estremo_fin;
                        r1 = punto_controllo - estremo_in;
                        r2 = punto_controllo - estremo_fin;

                        toll = 1e-10;
                        CP = cross(r1, r2);
                        CP_norm = dot(CP, CP);
                        if(CP_norm < toll)
                            CP_norm = toll;
                        end

                        U_ind = U_ind + (1/(4*pi))*dot(r0, (r1/norm(r1)) - (r2/norm(r2))).*CP./CP_norm;

                        % Induzione vortice portante
                        estremo_in = vortici{1,i_inducente}{j_inducente, k_inducente+1};
                        estremo_fin = vortici{1,i_inducente}{j_inducente, k_inducente};
                        r0 = estremo_in - estremo_fin;
                        r1 = punto_controllo - estremo_in;
                        r2 = punto_controllo - estremo_fin;

                        toll = 1e-10;
                        CP = cross(r1, r2);
                        CP_norm = dot(CP, CP);
                        if(CP_norm < toll)
                            CP_norm = toll;
                        end

                        U_ind = U_ind + (1/(4*pi))*dot(r0, (r1/norm(r1)) -(r2/norm(r2))).*CP./CP_norm;

                        % Induzione vortice semi-infinito dx
                        estremo_in = vortici_a_infinito{1,i_inducente}{j_inducente, k_inducente+1};
                        estremo_fin = vortici{1,i_inducente}{j_inducente, k_inducente+1};
                        r0 = estremo_in - estremo_fin;
                        r1 = punto_controllo - estremo_in;
                        r2 = punto_controllo - estremo_fin;

                        toll = 1e-10;
                        CP = cross(r1, r2);
                        CP_norm = dot(CP, CP);
                        if(CP_norm < toll)
                            CP_norm = toll;
                        end
                        U_ind = U_ind + (1/(4*pi))*dot(r0, (r1/norm(r1)) -(r2/norm(r2))).*CP./CP_norm;
                        
                        A(pannello_indotto, vortice_inducente) = dot(U_ind, normale);
                    end
                end
            end
        end
    end
end


% Costruzione termine noto
pannello_indotto = 0;

for i = 1:N_corpi
    for j_indotto = 1:discretizzazione_corda(i)
        for k_indotto = 1:2*discretizzazione_semiapertura_alare(i)
            pannello_indotto = pannello_indotto + 1;
            normale = normali_pannelli{1,i}{j_indotto, k_indotto};
            
            b(pannello_indotto) = -dot(U_inf, normale);
        end
    end
end


% Soluzione sistema lineare
sol = linsolve(A,b);
pannello = 0;

for i = 1:N_corpi
    gamma{i} = zeros(discretizzazione_corda(i), 2*discretizzazione_semiapertura_alare(i));

    for j = 1:discretizzazione_corda(i)
        for k = 1:2*discretizzazione_semiapertura_alare(i)
            pannello = pannello + 1;
            gamma{i}(j,k) = sol(pannello);
        end
    end
end
    
   



