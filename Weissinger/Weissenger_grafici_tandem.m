%% CODICE WEISSINGER CESSNA 172 GRAFICI COEFFICIENTI AERODINAMICI

clear
close all
clc

%% Dati velivolo
BSFC=7.126*1e-7; % lo voglio in [1/s]--> scegliere valore corretto
eta_p=0.9; % rendimento elica
Endurance_target=5*3600; % endurance voluta in [s]
W=10907.4; % peso velivolo in [N]


%% Dati iniziali

% velocità asintotica
U_inf = [51 0 0]';

% densità
rho = 0.9092; % Densità a 3000 metri

% angolo di sideslip (inserire in gradi)
beta = 0;
beta = deg2rad(beta);

% Rotazione della velocità
U_inf = [cos(beta) sin(beta) 0]' .* U_inf;

% Inserire numero di corpi
N_corpi = 2; % INSERIRE COME PRIMO CORPO QUELLO RISPETTO AL QUALE SI
% VOGLIONO ADIMENSIONALIZZARE I COEFFICIENTI AERODINAMICI TOTALI

corda_radice = [1.473 1.0311]';
diedro = [0 0]'; % inserire in gradi
diedro = deg2rad(diedro);

freccia = [0 0]'; % inserire in gradi
freccia = deg2rad(freccia);

rastremazione = [1 0.8]';
apertura_alare = [11 0.25*11]';

LE_posizione_X = [0 6]';
LE_posizione_Y = [0 0]';
LE_posizione_Z = [0.5 0]';

% angolo_di_rollio = [0];
alpha1_vect = [7.42, 0:1:10];
alpha2_vect = [-1.55, 0:-1:-5];
% angolo_di_beccheggio = [0];

% Discretizzazione pannellizzazione
discretizzazione_semiapertura_alare = [30 8]'; % numero panelli direzione apertura su SEMIALA
discretizzazione_corda = [10 5]';

% Inizializzazione sistema lineare
num_tot_pann = 0;
for i = 1:N_corpi
    num_tot_pann = num_tot_pann + 2*discretizzazione_corda(i)*discretizzazione_semiapertura_alare(i);
end

CL1=zeros(length(alpha1_vect),length(alpha2_vect));
CL2=zeros(length(alpha1_vect),length(alpha2_vect));
CD1=zeros(length(alpha1_vect),length(alpha2_vect));
CD2=zeros(length(alpha1_vect),length(alpha2_vect));
CL_tot_vect_1=zeros(length(alpha1_vect),length(alpha2_vect));
CL_trim_1=zeros(length(alpha1_vect),length(alpha2_vect));
CL_tot_vect_2=zeros(length(alpha1_vect),length(alpha2_vect));
CL_trim_2=zeros(length(alpha1_vect),length(alpha2_vect));


for a = 1:length(alpha1_vect)

    alpha1 = deg2rad(alpha1_vect(a));

    for z = 1:length(alpha2_vect)

        alpha2 = deg2rad(alpha2_vect(z));

        angolo_di_incidenza_aero(1,1) = alpha1;
        angolo_di_incidenza_aero(1,2) = alpha2;


        %% Creazione struttura ala
        estremi_ala = cell(1,N_corpi);

        mesh = cell(1,N_corpi);
        vortici = cell(1,N_corpi);
        punti_di_controllo = cell(1,N_corpi);
        punti_di_induzione = cell(1,N_corpi);
        normali_pannelli = cell(1,N_corpi);
        vortici_a_infinito = cell(1,N_corpi);
        superficie_alare = zeros(1,N_corpi);

        for i = 1:N_corpi

            %% CREAZIONE GEOMETRIA ALA
            % Parametri per costruzione ala

            semi_apertura_alare = apertura_alare(i) / 2;

            superficie_alare(i) = 2 * (semi_apertura_alare * corda_radice(i) * ( 1 + rastremazione(i) ) / 2);

            superficie_piana = superficie_alare(i) * cos(diedro(i));

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

            lunghezza_vortici_infinito = 50 * corda_radice; %% ATTENZIONE LUNGHEZZA VORTICI

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
        A = zeros(num_tot_pann, num_tot_pann);
        b = zeros(num_tot_pann, 1);

        gamma = cell(1, N_corpi);

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
            gamma{1,i} = zeros(discretizzazione_corda(i), 2*discretizzazione_semiapertura_alare(i));

            for j = 1:discretizzazione_corda(i)
                for k = 1:2*discretizzazione_semiapertura_alare(i)
                    pannello = pannello + 1;
                    gamma{1,i}(j,k) = sol(pannello);
                end
            end
        end


        %% Calcolo forze
        % Portanza
        gamma_sezione = cell(1,N_corpi);
        L_3D = cell(1,N_corpi);
        L_2D = cell(1,N_corpi);
        cL_3D = cell(1,N_corpi);
        delta_b_vect = cell(1,N_corpi);

        for i = 1:N_corpi
            L_3D_sezione = 0;
            for j = 1:2*discretizzazione_semiapertura_alare(i)
                gamma_sezione{1,i}{1,j} = sum(gamma{1,i}(:,j));
                gamma_i = cell2mat(gamma_sezione{1,i}(1,j));
                L_2D{1,i}{1,j} = rho*norm(U_inf)*gamma_i*cos(diedro(i));

                % deltab
                delta_b_vect{1,i}{1,j} = apertura_alare(i)/(2*discretizzazione_semiapertura_alare(i));
                delta_b = cell2mat(delta_b_vect{1,i}(1,j));

                % L_3D sezione
                L_2D_i = cell2mat(L_2D{1,i}(1,j));
                L_3D_sezione = L_3D_sezione + L_2D_i*delta_b;
            end

            % L_3D
            L_3D{1,i} = L_3D_sezione;
            cL_3D{1,i} = L_3D{1,i} / (0.5*rho*norm(U_inf)^2*superficie_alare(i));
        end


        % Resistenza
        punti_quarto_corda = cell(1, N_corpi);

        for i = 1 :N_corpi
            for k = 1:2*discretizzazione_semiapertura_alare(i)
                punto_medio_pannello_LE = (mesh{1,i}{1,k} - mesh{1,i}{1,k+1})/2 + mesh{1,i}{1,k+1};
                punto_medio_pannello_TE = (mesh{1,i}{length(discretizzazione_direzione_corda),k} - mesh{1,i}{length(discretizzazione_direzione_corda),k+1})/2 + mesh{1,i}{length(discretizzazione_direzione_corda),k+1};
                punto_quarto_corda = punto_medio_pannello_LE + 0.25*(punto_medio_pannello_TE-punto_medio_pannello_LE);
                punti_quarto_corda{1,i}{1,k} = punto_quarto_corda;
            end
        end

        num_profili_ala = cell(1, N_corpi);
        for i = 1:N_corpi
            num_profili_ala{1,i} = 2*discretizzazione_semiapertura_alare(i);
        end

        velocita_induzione_profili = cell(1, N_corpi);
        alpha_ind = cell(1, N_corpi);

        D_2D = cell(1,N_corpi);
        D_3D = cell(1,N_corpi);
        cD_3D = cell(1,N_corpi);

        for i = 1:N_corpi
            D_3D_sezione = 0;
            for k_indotto = 1:(num_profili_ala{1,i})

                v_ind = 0;
                punto_controllo = punti_quarto_corda{1,i}{1, k_indotto};
                normale= normali_pannelli{1,i}{1, k_indotto};


                for i_inducente = 1:N_corpi
                    for j_inducente = 1:discretizzazione_corda(i_inducente)
                        for k_inducente = 1:2*discretizzazione_semiapertura_alare(i_inducente)

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

                            gamma_ind = gamma{1,i_inducente}(j_inducente, k_inducente);
                            v_ind = v_ind + (gamma_ind/(4*pi))*dot(r0, (r1/norm(r1)) - (r2/norm(r2))).*CP./CP_norm;

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

                            gamma_ind = gamma{1,i_inducente}(j_inducente, k_inducente);
                            v_ind = v_ind + (gamma_ind/(4*pi))*dot(r0, (r1/norm(r1)) -(r2/norm(r2))).*CP./CP_norm;
                        end
                    end
                end
                velocita_induzione_profili{1,i}{1,k_indotto} = v_ind;

                % alpha indotto
                alpha_ind{1,i}{1, k_indotto} = atand(dot(v_ind, normale)/norm(U_inf));
                alpha_i = alpha_ind{1,i}{1, k_indotto};

                % drag 2D
                D_2D{1,i}{1,k_indotto} = abs(L_2D{1,i}{1, k_indotto})*sind(abs(alpha_i));
                D_2D_i = D_2D{1,i}{1,k_indotto};

                % deltab
                delta_b = cell2mat(delta_b_vect{1,i}(1,k_indotto));

                % drag 3D
                D_3D_sezione = D_3D_sezione + delta_b*D_2D_i;
                D_3D{1,i} = D_3D_sezione;
                cD_3D{1,i} = D_3D_sezione / (0.5*rho*norm(U_inf)^2*superficie_alare(i));
            end
        end

        %% Calcolo Cl totale e Cd totale
        L_tot = sum(cell2mat(L_3D(1,:)));
        CL_tot=L_tot/(0.5*rho*(norm(U_inf)^2*superficie_alare(1)));
        D_tot = sum(cell2mat(D_3D(1,:)));
        CD_tot=D_tot/(0.5*rho*(norm(U_inf)^2*superficie_alare(1)));


        CL1(a,z)=cL_3D{1, 1};
        CL2(a,z)=cL_3D{1, 2};
        CD1(a,z)=cD_3D{1, 1};
        CD2(a,z)=cD_3D{1, 2};

        CL_tot_vect_1(a,z)=CL_tot;
        CL_trim_1(a,z)=2*W/(rho*norm(U_inf)^2*superficie_alare(1));
        CL_tot_vect_2(a,z)=CL_tot;
        CL_trim_2(a,z)=2*W/(rho*norm(U_inf)^2*superficie_alare(1));
    end

end

%% Grafico con alpha1 variabile e alpha2 fisso ad esempio
% Rappresento, in pratica, una colonna di Cl corrispondente all'alpha2
% voluto (esempio: prendo la prima colonna, cioè alpha2=2)--> ho il grafico
% di Cl e Cd in funzione di ALPHA1_vect, con ALPHA2 fissato
figure
plot(alpha1_vect,CL1(:,1)','LineWidth',5)
hold on
plot(alpha1_vect,CL2(:,1)','LineWidth',5)
hold on
plot(alpha1_vect,CD1(:,1)','LineWidth',5)
hold on
plot(alpha1_vect,CD2(:,1)','LineWidth',5)
legend('C_{L_1}','C_{L_2}','C_{D_1}','C_{D_2}')
grid on
xlabel('\alpha_1', 'FontSize', 40)
ylabel('Coeff aerodinamici', 'FontSize', 40)
ax = gca;
ax.FontSize = 40;
title('Coefficienti aerodinamici con \alpha_2 fissato', 'FontSize', 50)

%% Grafico con alpha2 variabile e alpha1 fisso ad esempio
% Rappresento, in pratica, una riga di Cl corrispondente all'alpha1
% voluto (esempio: prendo la prima riga, cioè alpha1=4)--> ho il grafico
% di Cl e Cd in funzione di ALPHA2_vect, con ALPHA1 fissato
figure
plot(alpha2_vect,CL1(1,:)','LineWidth',5)
hold on
plot(alpha2_vect,CL2(1,:)','LineWidth',5)
hold on
plot(alpha2_vect,CD1(1,:)','LineWidth',5)
hold on
plot(alpha2_vect,CD2(1,:)','LineWidth',5)
legend('C_{L_1}','C_{L_2}','C_{D_1}','C_{D_2}')
grid on
xlabel('\alpha_2', 'FontSize', 40)
ylabel('Coeff aerodinamici', 'FontSize', 40)
ax = gca;
ax.FontSize = 40;
title('Coefficienti aerodinamici con \alpha_1 fissato', 'FontSize', 50)


%% Grafico alpha1 trim
figure
plot(alpha1_vect, CL_trim_1(:,1), 'lineStyle', '--', 'LineWidth',5)
hold on
plot(alpha1_vect, CL_tot_vect_1(:,1), 'LineWidth',5)
xline(7.42, 'lineStyle', '-.','LineWidth',5)
plot(7.42, CL_trim_1(1,1), '.', 'MarkerSize', 50, 'Color', 'k')
legend('C_{L_{trim}}','C_{L_{tot}}', '\alpha_{1_{trim}}')
xlabel('\alpha_1', 'FontSize', 40)
ylabel('C_{L_{tot}}', 'FontSize', 40)
ax = gca;
ax.FontSize = 40;
title('\alpha_1 di trim', 'FontSize', 50)
grid on


%% Grafico alpha2 trim
figure
plot(alpha2_vect, CL_trim_2(1,:), 'LineWidth',5)
hold on
plot(alpha2_vect, CL_tot_vect_2(1,:), 'LineWidth',5)
xline(-1.65, 'lineStyle', '-.','LineWidth',5)
plot(-1.65, CL_trim_1(1,1), '.', 'MarkerSize', 50, 'Color', 'k')
xlabel('\alpha_2', 'FontSize', 40)
ylabel('C_{L_{tot}}', 'FontSize', 40)
legend('C_{L_{trim}}','C_{L_{tot}}', '\alpha_{2_{trim}}', 'FontSize', 50)
ax = gca;
ax.FontSize = 40;
title('\alpha_2 di trim', 'FontSize', 50)
grid on

