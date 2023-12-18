% Il risultato dell'effetto ground è un aumento della pendenza della curva
% Cl-alpha; ciò risulta in una diminuzione del Cl_ZL (zero-lift), ma un
% aumento del Cl per alpha>1-2°


%% Pulizia
clear
clc
close all
addpath mat_functions

%% Input
% Osservo che, poichè il NACA 0012 è simmetrico, se alpha=0° non riesco a
% generare portanza--> ho circolazione complessiva nulla sul profilo (ciò
% risulta nel fatto che l'intensità dei vortici x_v sia nulla).

h=0.5; % distanza tra profilo e suolo (semidistanza verticale tra i 2 
% profili da utilizzare con metodo delle immagini

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = -10;   % opposto dell'inclinazione profilo rispetto alla corr uniforme (incidenza)

% La corrente uniforme deve rimanere orizzontale, perchè sennò non potrei
% creare correttamente lo specchio (non riuscirei a creare una parete
% orizzontale con una U_inf inclinata)
U_inf_x = U_inf;
U_inf_y = 0;

U_inf = [U_inf_x; U_inf_y]; % vettore di U_inf in sist globale (X,Y)
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal); % versore normale ad U_inf

%% Creazione profilo
TestCase = 0;
NCorpi = 1;  % Numero di corpi da analizzare

CodiceProfilo = cell(NCorpi, 1);
CodiceProfilo{1} = '23012';
Chord = [1];
NPannelli = [101];

i=1;
[x,y]=createProfile(CodiceProfilo{i},NPannelli(i),Chord(i), 1);

x=x.*cosd(alpha)-y.*sind(alpha);
y=x.*sind(alpha)+y.*cosd(alpha);

naca0012nodes_profilo=[x y+h];
naca0012nodes_specchio=[x -h-y];
%% Plot profilo
plot(naca0012nodes_profilo(:,1),naca0012nodes_profilo(:,2),'Linewidth',2)
axis equal
hold on
plot(naca0012nodes_profilo(:,1),naca0012nodes_profilo(:,2),'.','MarkerSize',10,'Color','r')
hold on
plot(naca0012nodes_specchio(:,1),naca0012nodes_specchio(:,2),'Linewidth',2)
hold on
plot(naca0012nodes_specchio(:,1),naca0012nodes_specchio(:,2),'.','MarkerSize',10,'Color','r')
grid on

%% Creo pannelli profilo
% Numero di pannelli
N_pann=size(naca0012nodes_profilo,1)-1;
% Vettori dei pannelli (uso senso antiorario a partire dal nodo 1 della
% lista)--> li salvo come tanti vettori riga
pan_vect=zeros(N_pann,2);
centre_vect=zeros(N_pann,2);
estr_1_vect=zeros(N_pann,2);
estr_2_vect=zeros(N_pann,2);
tangent_versor=zeros(N_pann,2);
normal_versor=zeros(N_pann,2);
lungh_vect=zeros(N_pann,1);
for i=1:N_pann
    estr_1_vect(i,:)=naca0012nodes_profilo(i,:);
    estr_2_vect(i,:)=naca0012nodes_profilo(i+1,:);
    pan_vect(i,:)=naca0012nodes_profilo(i+1,:)-naca0012nodes_profilo(i,:);
    centre_vect(i,:)=naca0012nodes_profilo(i,:)+pan_vect(i,:)/2;
    tangent_versor(i,:)=pan_vect(i,:)/norm(pan_vect(i,:));
    normal_versor(i,1)=-tangent_versor(i,2);
    normal_versor(i,2)=tangent_versor(i,1);
    lungh_vect(i,1)=norm(pan_vect(i,:));
end
%% Creo pannelli specchio
% Numero di pannelli
N_pann_s=size(naca0012nodes_specchio,1)-1;
% Vettori dei pannelli (uso senso antiorario a partire dal nodo 1 della
% lista)--> li salvo come tanti vettori riga
pan_vect_s=zeros(N_pann_s,2);
centre_vect_s=zeros(N_pann_s,2);
estr_1_vect_s=zeros(N_pann_s,2);
estr_2_vect_s=zeros(N_pann_s,2);
tangent_versor_s=zeros(N_pann_s,2);
normal_versor_s=zeros(N_pann_s,2);
lungh_vect_s=zeros(N_pann_s,1);
for i=1:N_pann_s
    estr_1_vect_s(i,:)=naca0012nodes_specchio(i,:);
    estr_2_vect_s(i,:)=naca0012nodes_specchio(i+1,:);
    pan_vect_s(i,:)=naca0012nodes_specchio(i+1,:)-naca0012nodes_specchio(i,:);
    centre_vect_s(i,:)=naca0012nodes_specchio(i,:)+pan_vect_s(i,:)/2;
    tangent_versor_s(i,:)=pan_vect_s(i,:)/norm(pan_vect_s(i,:));
    normal_versor_s(i,1)=-tangent_versor_s(i,2);
    normal_versor_s(i,2)=tangent_versor_s(i,1);
    lungh_vect_s(i,1)=norm(pan_vect_s(i,:));
end
%% Rappresento i centri dei pannelli
plot(naca0012nodes_profilo(:,1),naca0012nodes_profilo(:,2),'Linewidth',2)
axis equal
grid on
hold on
plot(naca0012nodes_profilo(:,1),naca0012nodes_profilo(:,2),'.','MarkerSize',7,'Color','r')
hold on
plot(centre_vect(:,1),centre_vect(:,2),'.','MarkerSize',7,'Color','g')

%% Matrice di rotazione profilo (sist locale--> sist globale)
Loc_to_Glob_matrix=zeros(2*N_pann,2);
for j=1:N_pann
    Loc_to_Glob_matrix((2*j-1:2*j),1)=tangent_versor(j,:)';
    Loc_to_Glob_matrix((2*j-1:2*j),2)=normal_versor(j,:)';
end
% La matrice del pannello i-esimo sarà: R_i=Loc_to_Glob_matrix((i,i+1),:)

%% Matrice di rotazione specchio (sist locale--> sist globale)
Loc_to_Glob_matrix_s=zeros(2*N_pann_s,2);
for j=1:N_pann_s
    Loc_to_Glob_matrix_s((2*j-1:2*j),1)=tangent_versor_s(j,:)';
    Loc_to_Glob_matrix_s((2*j-1:2*j),2)=normal_versor_s(j,:)';
end
%% Scrivo parte quadrata matrice A del sistema lineare
% Velocità che il profilo induce su se stesso
a11 = zeros(N_pann,1);
for i = 1: N_pann
    for j = 1: N_pann
        R_j=Loc_to_Glob_matrix((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect(i,:)',estr_1_vect(j,:)',...
            estr_2_vect(j,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect(i,:)',estr_1_vect(j,:)',...
            estr_2_vect(j,:)',R_j,R_j');  
        A11(i,j) = Us_j_ri(1,1)*normal_versor(i,1) + Us_j_ri(2,1)*normal_versor(i,2);
        a11(i,1) = a11(i,1) + Uv_j_ri(1,1)*normal_versor(i,1) + Uv_j_ri(2,1)*normal_versor(i,2);
    end
end
% Velocità che lo specchio induce sul profilo
for i = 1: N_pann_s
    for j = 1: N_pann_s
        R_j_s=Loc_to_Glob_matrix_s((2*j-1:2*j),:);
        [Us_j_ri_s] = ViSorgente(centre_vect(i,:)',estr_1_vect_s(j,:)',...
            estr_2_vect_s(j,:)',R_j_s,R_j_s');
        [Uv_j_ri_s] = ViVortice(centre_vect(i,:)',estr_1_vect_s(j,:)',...
            estr_2_vect_s(j,:)',R_j_s,R_j_s');
        A11(i,j) = A11(i,j) + (Us_j_ri_s(1,1)*normal_versor(i,1) + Us_j_ri_s(2,1)*normal_versor(i,2));
        a11(i,1) = a11(i,1) - (Uv_j_ri_s(1,1)*normal_versor(i,1) + Uv_j_ri_s(2,1)*normal_versor(i,2));
        % qui il segno è "meno" perchè il vortice specchiato cambia verso
    end
end

%% Completo la scrittura di A (usando condizione di Kutta)
% Induzione profilo su se stesso (nel primo e nell'ultimo punto di
% controllo)
c11=0;
for j = 1:N_pann
    R_j=Loc_to_Glob_matrix((2*j-1:2*j),:);
    [Us_j_r1] = ViSorgente(centre_vect(1,:)',estr_1_vect(j,:)',...
        estr_2_vect(j,:)',R_j,R_j');
    [Uv_j_r1] = ViVortice(centre_vect(1,:)',estr_1_vect(j,:)',...
        estr_2_vect(j,:)',R_j,R_j');
    C11(1,j) = dot(Us_j_r1,tangent_versor(1,:)');
    c11 = c11 + dot(Uv_j_r1,tangent_versor(1,:)');
    [Us_j_rN] = ViSorgente(centre_vect(N_pann,:)',estr_1_vect(j,:)',...
        estr_2_vect(j,:)',R_j,R_j');
    [Uv_j_rN] = ViVortice(centre_vect(N_pann,:)',estr_1_vect(j,:)',...
        estr_2_vect(j,:)',R_j,R_j');
    C11(1,j) = C11(1,j) + dot(Us_j_rN,tangent_versor(N_pann,:)'); 
    c11 = c11 + dot(Uv_j_rN,tangent_versor(N_pann,:)');    
end
% Induzione specchio su profilo (nel primo e nell'ultimo punto di
% controllo)
for j = 1:N_pann_s
    R_j_s=Loc_to_Glob_matrix_s((2*j-1:2*j),:);
    [Us_j_r1] = ViSorgente(centre_vect(1,:)',estr_1_vect_s(j,:)',...
        estr_2_vect_s(j,:)',R_j_s,R_j_s');
    [Uv_j_r1] = ViVortice(centre_vect(1,:)',estr_1_vect_s(j,:)',...
        estr_2_vect_s(j,:)',R_j_s,R_j_s');
    C11(1,j) = C11(1,j)+dot(Us_j_r1,tangent_versor(1,:)');
    c11 = c11 - dot(Uv_j_r1,tangent_versor(1,:)');
    [Us_j_rN] = ViSorgente(centre_vect(N_pann,:)',estr_1_vect_s(j,:)',...
        estr_2_vect_s(j,:)',R_j_s,R_j_s');
    [Uv_j_rN] = ViVortice(centre_vect(N_pann,:)',estr_1_vect_s(j,:)',...
        estr_2_vect_s(j,:)',R_j_s,R_j_s');
    C11(1,j) = C11(1,j) + dot(Us_j_rN,tangent_versor(N_pann,:)'); 
    c11 = c11 - dot(Uv_j_rN,tangent_versor(N_pann,:)');  
end

%% Assemblo matrice A
A=zeros(N_pann+1);
A(1:N_pann,1:N_pann)=A11;
A(1:N_pann,N_pann+1)=a11;
A(N_pann+1,1:N_pann)=C11;
A(N_pann+1,N_pann+1)=c11;

%% Scrivo il termine noto del sistema lineare
b=zeros(N_pann+1,1);
for i=1:N_pann
    b(i)=-dot(U_inf',normal_versor(i,:));
end
b(end,1)=-dot(U_inf',tangent_versor(1,:))-dot(U_inf',tangent_versor(end,:));

%% Risolvo il sistema lineare: trovo x_vect
x_vect=linsolve(A,b);

%% Calcolo velocità complessiva nei vari punti di controllo
% Poichè il Cp mi interessa sulla superficie del profilo, lo calcolo con la
% velocità che ho sulla superficie dei pannelli (in particolare, quella
% nei punti di controllo)--> sulla superficie dei pannelli ho solo velocità
% tangenziale (da condizione di impenetrabilità)

u_ri=zeros(N_pann,1); % le righe sono la compon x di U_ri 
v_ri=zeros(N_pann,1);
% airfoil 1
for i =1:N_pann
    u_ri(i)=U_inf(1); 
    v_ri(i)=U_inf(2);
    for j = 1:N_pann
        R_j=Loc_to_Glob_matrix((2*j-1:2*j),:);
        % velocità indotta dal profilo su se stesso
        [Us_j_ri] = ViSorgente(centre_vect(i,:)',estr_1_vect(j,:)',estr_2_vect(j,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect(i,:)',estr_1_vect(j,:)',estr_2_vect(j,:)',R_j,R_j');
        u_ri(i) = u_ri(i) + x_vect(j)*Us_j_ri(1,1) + x_vect(end)*Uv_j_ri(1,1);
        v_ri(i) = v_ri(i) + x_vect(j)*Us_j_ri(2,1) + x_vect(end)*Uv_j_ri(2,1);
        % Velocità indotta dallo specchio sul pannello
        R_j_s=Loc_to_Glob_matrix_s((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect(i,:)',estr_1_vect_s(j,:)',estr_2_vect_s(j,:)',R_j_s,R_j_s');
        [Uv_j_ri] = ViVortice(centre_vect(i,:)',estr_1_vect_s(j,:)',estr_2_vect_s(j,:)',R_j_s,R_j_s');
        u_ri(i) = u_ri(i) - (-x_vect(j)*Us_j_ri(1,1)+x_vect(end)*Uv_j_ri(1,1));
        v_ri(i) = v_ri(i) - (-x_vect(j)*Us_j_ri(2,1)+x_vect(end)*Uv_j_ri(2,1));        
    end
end

% Calcolo la velocità nei vari punti di controllo
U_ri_vect=[u_ri v_ri];

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

%% Calcolo circolazione totale profilo e CL (uso Kutta-Joukowsky)
Gamma=sum(x_vect(end).*lungh_vect);
% Ottengo il Cl adimensionalizzando con 0.5*rho*U_inf^2*c la portanza per
% unità di apertura ottenuta con Kutta-Joukowsky (cioè -rho*U_inf*Gamma)
Cl=2*Gamma/norm(U_inf);
% Oppure moltiplico il Cp per le varie lunghezze dei pannelli (come se
% integrassi la pressione sui pannelli) e poi adimensionalizzo con la corda
% "c" del profilo (ricordo che c=1)
Cl_p=0;
for i=1:N_pann
    Cl_p=Cl_p-dot(((1/Chord)*Cp(i,1)*lungh_vect(i).*normal_versor(i,:)),U_inf_normal);
end
%% Calcolo momento di beccheggio rispetto al polo "bordo d'attacco"
Cm=0;
for i=1:N_pann
    rc=[centre_vect(i,1) centre_vect(i,2)  0];
    nc=[normal_versor(i,1) normal_versor(i,2) 0];
    Cm=Cm-((1/Chord^2)*Cp(i,1)*lungh_vect(i)).*dot((cross(rc,nc)),[0 0 -1]);
    % il Cm è positivo a cabrare
end












