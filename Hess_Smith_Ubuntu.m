clear
close all
clc

%% Input
% Osservo che, poichè il NACA 0012 è simmetrico, se alpha=0° non riesco a
% generare portanza--> ho circolazione complessiva nulla sul profilo (ciò
% risulta nel fatto che l'intensità dei vortici x_v sia nulla).

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 0;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha)); % U_inf proiettata su X (sist globale)
U_inf_y = U_inf * sin(deg2rad(alpha)); % U_inf proiettata su Y (sist globale)

U_inf = [U_inf_x; U_inf_y]; % vettore di U_inf in sist globale (X,Y)
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal); % versore normale ad U_inf

TestCase = 0;
NCorpi = 1;  % Numero di corpi da analizzare

CodiceProfilo = cell(NCorpi, 1);
CodiceProfilo{1} = '0012';
Chord = [1];
NPannelli = [102];

LE_X_Position = [0];
LE_Y_Position = [0];
% Impongo che il bordo d'attacco sia nell'origine del sist di rif globale
% (X,Y)--> ricordo che il profilo su XFOIL è definito per X che va da 0 a 1.
LE_X_Position = [0];
LE_Y_Position = [0];

%% Creazione profilo

% numero profilo:
i=1;

[x,y]=createProfile(CodiceProfilo{i},NPannelli(i),Chord(i));

Corpi{i}.x=x;
Corpi{i}.y=y;

naca_nodes = [x, y];

% Devo togliere il gap sul bordo d'uscita usato da XFoil
naca_nodes(1,:)=(naca_nodes(1,:)+naca_nodes(end,:))/2;
naca_nodes(end,:)=naca_nodes(1,:);

%% Plot profilo
plot(naca_nodes(:,1),naca_nodes(:,2),'Linewidth',2)
grid on
axis equal
hold on
plot(naca_nodes(:,1),naca_nodes(:,2),'.','MarkerSize',10,'Color','r')

%% Creo i pannelli e trovo le grandezze geometriche importanti
% Numero di pannelli
N_pann=size(naca_nodes,1)-1;
% Vettori dei pannelli (uso senso antiorario a partire dal nodo 1 della
% lista)--> li salvo come tanti vettori riga
pan_vect=zeros(N_pann,2);
centre_vect=zeros(N_pann,2);
estr_1_vect=zeros(N_pann,2);
estr_2_vect=zeros(N_pann,2);
tangent_versor=zeros(N_pann,2);
normal_versor=zeros(N_pann,2);
lungh_vect=zeros(N_pann,2);
for i=1:N_pann
    estr_1_vect(i,:)=naca_nodes(i,:);
    estr_2_vect(i,:)=naca_nodes(i+1,:);
    pan_vect(i,:)=naca_nodes(i+1,:)-naca_nodes(i,:);
    centre_vect(i,:)=naca_nodes(i,:)+pan_vect(i,:)/2;
    tangent_versor(i,:)=pan_vect(i,:)/norm(pan_vect(i,:));
    normal_versor(i,1)=-tangent_versor(i,2);
    normal_versor(i,2)=tangent_versor(i,1);
    lungh_vect(i,:)=norm(pan_vect(i,:));
end

%% Rappresento i centri dei pannelli
plot(naca_nodes(:,1),naca_nodes(:,2),'Linewidth',2)
axis equal
grid on
hold on
plot(naca_nodes(:,1),naca_nodes(:,2),'.','MarkerSize',7,'Color','r')
hold on
plot(centre_vect(:,1),centre_vect(:,2),'.','MarkerSize',7,'Color','g')

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

%% Calcolo Cp nei vari punti di controllo
Cp=zeros(N_pann,1);
for i=1:N_pann
    Cp(i,1)=1-norm(U_ri_vect(i,:))^2/norm(U_inf)^2;
end
% Il vettore Cp è calcolato con l'ordine orario dei punti di controllo, 
% a partire dal punto di controllo del pannello 1

figure
plot(centre_vect(:, 1), -Cp(:, 1));
axis equal
grid on
grid minor
xlabel('x');
ylabel('-Cp');
