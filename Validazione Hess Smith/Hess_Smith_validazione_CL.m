%% Pulizia
clear
clc
close all
addpath mat_functions

%% Input
% Osservo che, poichè il NACA 0012 è simmetrico, se alpha=0° non riesco a
% generare portanza--> ho circolazione complessiva nulla sul profilo (ciò
% risulta nel fatto che l'intensità dei vortici x_v sia nulla).

alpha_vect=0:0.5:6;
Cl_vect=[];

for l=1:length(alpha_vect)
    alpha=alpha_vect(l);
    U_inf = 1;  % Velocità all'infinito [m/s]
    U_inf_x = U_inf * cos(deg2rad(alpha)); % U_inf proiettata su X (sist globale)
    U_inf_y = U_inf * sin(deg2rad(alpha)); % U_inf proiettata su Y (sist globale)
    rho=1.225;
    
    U_inf = [U_inf_x; U_inf_y]; % vettore di U_inf in sist globale (X,Y)
    U_inf_normal = [-U_inf(2); U_inf(1)];
    U_inf_normal = U_inf_normal ./ norm(U_inf_normal); % versore normale ad U_inf
    
    % Impongo che il bordo d'attacco sia nell'origine del sist di rif globale
    % (X,Y)--> ricordo che il profilo su XFOIL è definito per X che va da 0 a 1.
    LE_X_Position = [0];
    LE_Y_Position = [0];
    
    %% Creazione profilo
    TestCase = 0;
    NCorpi = 1;  % Numero di corpi da analizzare
    
    CodiceProfilo = cell(NCorpi, 1);
    CodiceProfilo{1} = '0012';
    Chord = [1];
    NPannelli = [101];
    
    i=1;
    [x,y]=createProfile(CodiceProfilo{i},NPannelli(i),Chord(i));
    
    naca0012nodes=[x y];
    %% Plot profilo
    % plot(naca0012nodes(:,1),naca0012nodes(:,2),'Linewidth',2)
    % axis equal
    % hold on
    % plot(naca0012nodes(:,1),naca0012nodes(:,2),'.','MarkerSize',10,'Color','r')
    % hold off
    
    %% Creo i pannelli e trovo le grandezze geometriche importanti
    % Numero di pannelli
    N_pann=size(naca0012nodes,1)-1;
    % Vettori dei pannelli (uso senso orario a partire dal nodo 1 della
    % lista)--> li salvo come tanti vettori riga
    pan_vect=zeros(N_pann,2);
    centre_vect=zeros(N_pann,2);
    estr_1_vect=zeros(N_pann,2);
    estr_2_vect=zeros(N_pann,2);
    tangent_versor=zeros(N_pann,2);
    normal_versor=zeros(N_pann,2);
    lungh_vect=zeros(N_pann,1);
    for i=1:N_pann
        estr_1_vect(i,:)=naca0012nodes(i,:);
        estr_2_vect(i,:)=naca0012nodes(i+1,:);
        pan_vect(i,:)=naca0012nodes(i+1,:)-naca0012nodes(i,:);
        centre_vect(i,:)=naca0012nodes(i,:)+pan_vect(i,:)/2;
        tangent_versor(i,:)=pan_vect(i,:)/norm(pan_vect(i,:));
        normal_versor(i,1)=-tangent_versor(i,2);
        normal_versor(i,2)=tangent_versor(i,1);
        lungh_vect(i,1)=norm(pan_vect(i,:));
    end
    %% Rappresento i centri dei pannelli
    % plot(naca0012nodes(:,1),naca0012nodes(:,2),'Linewidth',2)
    % axis equal
    % grid on
    % hold on
    % plot(naca0012nodes(:,1),naca0012nodes(:,2),'.','MarkerSize',7,'Color','r')
    % hold on
    % plot(centre_vect(:,1),centre_vect(:,2),'.','MarkerSize',7,'Color','g')
    
    %% Matrice di rotazione da sist di rif locale a globale: U=R*u.
    Loc_to_Glob_matrix=zeros(2*N_pann,2);
    for j=1:N_pann
        Loc_to_Glob_matrix((2*j-1:2*j),1)=tangent_versor(j,:)';
        Loc_to_Glob_matrix((2*j-1:2*j),2)=normal_versor(j,:)';
    end
    % La matrice del pannello i-esimo sarà: R_i=Loc_to_Glob_matrix((i,i+1),:)
    
    %% Velocità data da sorgenti in tutti i centri dei pannelli (sist globale)
    Us_j_ri_matr=zeros(N_pann,2*N_pann);
    for j=1:N_pann
        Estremo_1_j=estr_1_vect(j,:)';
        Estremo_2_j=estr_2_vect(j,:)';
        R_j=Loc_to_Glob_matrix(((2*j-1):2*j),:);
        for i=1:N_pann
            Centro_i=centre_vect(i,:)';
            [Us_j_ri]=ViSorgente(Centro_i, Estremo_1_j, Estremo_2_j, R_j, R_j');
            Us_j_ri_matr(i,(2*j-1:2*j))=Us_j_ri'; % li salvo come tanti vettori riga
        end
    end
    % La velocità Us_j(ri) è il termine Us_j_ri_matr(i,(2*j-1:2*j))--> la riga
    % indica il punto in cui induco la velocità, mentre la colonna indica il
    % pannello che genera tale velocità tramite le proprie sorgenti.
    
    %% Velocità data da vortici in tutti i centri dei pannelli (sist globale)
    Uv_j_ri_matr=zeros(N_pann,2*N_pann);
    for j=1:N_pann
        Estremo_1_j=estr_1_vect(j,:)';
        Estremo_2_j=estr_2_vect(j,:)';
        R_j=Loc_to_Glob_matrix((2*j-1:2*j),:);
        for i=1:N_pann
            Centro_i=centre_vect(i,:)';
            [Uv_j_ri]=ViVortice(Centro_i, Estremo_1_j, Estremo_2_j, R_j, R_j');
            Uv_j_ri_matr(i,(2*j-1:2*j))=Uv_j_ri'; % li salvo come tanti vettori riga
        end
    end
    % La velocità Uv_j(ri) è il termine Uv_j_ri_matr(i,(2*j-1:2*j))--> la riga
    % indica il punto in cui induco la velocità, mentre la colonna indica il
    % pannello che genera tale velocità tramite i propri vortici.
    
    %% Scrivo il termine noto del sistema lineare
    b=zeros(N_pann+1,1);
    for i=1:N_pann
        b(i)=-dot(U_inf',normal_versor(i,:));
    end
    b(end,1)=-dot(U_inf',tangent_versor(1,:))-dot(U_inf',tangent_versor(end,:));
    
    %% Scrivo la matrice A del sistema lineare
    A=zeros(N_pann+1,N_pann+1);
    % Calcolo A(i,j)
    for i=1:N_pann
        for j=1:N_pann
            A(i,j)=dot(Us_j_ri_matr(i,(2*j-1:2*j)),normal_versor(i,:));
        end
    end
    
    % Calcolo A(N+1,j)
    for j=1:N_pann
        A(end,j)=dot(Us_j_ri_matr(1,(2*j-1:2*j)),tangent_versor(1,:))+...
        dot(Us_j_ri_matr(end,(2*j-1:2*j)),tangent_versor(end,:));
    end
    
    % Calcolo A(N+1,N+1)
    for j=1:N_pann
        A(end,end)=A(end,end)+dot(Uv_j_ri_matr(1,(2*j-1:2*j)),tangent_versor(1,:))+...
               dot(Uv_j_ri_matr(end,(2*j-1:2*j)),tangent_versor(end,:));
    end
    
    % Calcolo A(i,N+1)
    for i=1:N_pann
        for j=1:N_pann
            A(i,end)=A(i,end)+dot(Uv_j_ri_matr(i,(2*j-1:2*j)),normal_versor(i,:));
        end
    end
    
    %% Risolvo il sistema lineare: trovo x_vect
    x_vect=linsolve(A,b);
    
    %% Calcolo velocità complessiva nei vari punti di controllo
    % Poichè il Cp mi interessa sulla superficie del profilo, lo calcolo con la
    % velocità che ho sulla superficie dei pannelli (in particolare, quella
    % nei punti di controllo)--> sulla superficie dei pannelli ho solo velocità
    % tangenziale (da condizione di impenetrabilità)
    
    % Calcolo la velocità nei vari punti di controllo
    U_ri_vect=zeros(N_pann,2);
    for i=1:N_pann
        for j=1:N_pann
            U_ri_vect(i,:)=U_ri_vect(i,:)+Us_j_ri_matr(i,(2*j-1:2*j))*x_vect(j)+...
                 Uv_j_ri_matr(i,(2*j-1:2*j))*x_vect(end);
        end
    end
    U_ri_vect(:,1)=U_inf(1,1)+U_ri_vect(:,1);
    U_ri_vect(:,2)=U_inf(2,1)+U_ri_vect(:,2);
    
    % Le righe di U_ri_vect rappresentano la velocità globale nel punto di
    % controllo del pannello i-esimo.
    
    %% Controllo se la velocità globale è tangente ad ogni punto di controllo
    prod_scal=[];
    for i=1:N_pann
        prod_scal(i)=dot(U_ri_vect(i,:),normal_versor(i,:));
    end
    % ok
    %% Controllo condizione di Kutta
    zero=dot(U_ri_vect(1,:),tangent_versor(1,:))+dot(U_ri_vect(end,:),...
        tangent_versor(end,:));
    % ok
    %% Calcolo Cp nei vari punti di controllo
    Cp=zeros(N_pann,1);
    for i=1:N_pann
        Cp(i,1)=1-norm(U_ri_vect(i,:))^2/norm(U_inf)^2;
    end
    % Il vettore Cp è calcolato con l'ordine orario dei punti di controllo, 
    % a partire dal punto di controllo del pannello 1
    
    %% Plot del Cp
    % figure
    % plot(centre_vect((1:(N_pann+1)/2),1),-Cp((1:(N_pann+1)/2),1),'r','LineWidth',2)
    % hold on
    % plot(centre_vect((N_pann+1)/2:end,1),-Cp((N_pann+1)/2:end,1),'b','LineWidth',2)
    % legend('-Cp ventre profilo','-Cp dorso profilo')
    % grid on
    % xlim([0 1])
    % hold off
    
    %% Calcolo circolazione totale profilo e CL (uso Kutta-Joukowsky)
    Gamma=sum(x_vect(end).*lungh_vect);
    Cl=2*Gamma/norm(U_inf);
    % Oppure moltiplico il Cp per le varie lunghezze dei pannelli (come se
    % integrassi la pressione sui pannelli) e poi adimensionalizzo con la corda
    % "c" del profilo (ricordo che c=1)
    Cl_p=0;
    for i=1:N_pann
        Cl_p=Cl_p-dot((1/Chord)*Cp(i,1)*lungh_vect(i).*normal_versor(i,:),U_inf_normal);
    end
    
    %% Calcolo momento di beccheggio rispetto al polo "(0.25,0)"
    Cm=0;
    for i=1:N_pann
        rc=[centre_vect(i,1)-0.25 centre_vect(i,2)  0];
        nc=[normal_versor(i,1) normal_versor(i,2) 0];
        Cm=Cm-((1/Chord^2)*Cp(i,1)*lungh_vect(i)).*dot((cross(rc,nc)),[0 0 -1]);
        % il Cm è positivo a cabrare
    end
    Cl_vect(end+1)=Cl_p;
end

%% Plot Cl-alpha naca di XFoil e del profilo
figure
polareNaca0012 = importfile('polareNaca0012');
plot(polareNaca0012(1:13,1),polareNaca0012(1:13,2),'LineWidth',5,'Color','r','LineStyle','-')
hold on
plot(alpha_vect,Cl_vect,'Linewidth',5,'Color','g')
title('Confronto tra i grafici del Cl(\alpha)','FontSize',50)
grid on
legend('XFoil','Hess Smith','Fontsize',40)
ax=gca;
ax.FontSize=40;
xlim([0 6])
hold off
xlabel('alpha','FontSize',40)
ylabel('C_L','FontSize',40)

%% Errore relativo massimo tra il Cl di Hess Smith e quello di XFoil
err_rel_max_Cl=max(abs(Cl_vect-polareNaca0012(1:13,2)')./polareNaca0012(1:13,2)')*100;
% l'errore relativo massimo per alpha tra 0 e 6 è di 1.935%




