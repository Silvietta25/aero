clear all
close all
addpath mat_functions

%% Distanza tra i profili in tandem
x12 = 1.5;  % distance between forward airfoil LE and back airfoil LE
y12 = 0;

%% Input
U = 1;
alpha1 = -5*pi/180; % è opposta all'inclinazione del profilo 1 rispetto ad asse X
alpha2 =-10*pi/180; % è opposta all'inclinazione del profilo 2 rispetto ad asse X

U_inf(1)=U; 
U_inf(2)=0;

%% Dati profili
TestCase = 0;
NCorpi = 2;  % Numero di corpi da analizzare
CodiceProfilo = cell(NCorpi, 1);
CodiceProfilo{1} = '0012';
CodiceProfilo{2} = '0012';
Chord = [1 1];
N_pann1 = 101;
N_pann2 = 101;

%% Creazione profilo 1
i=1;
[x_nodes1,y_nodes1]=createProfile(CodiceProfilo{i},N_pann1,Chord(i));
x_nodes1=x_nodes1.*cos(alpha1)-y_nodes1.*sin(alpha1);
y_nodes1=x_nodes1.*sin(alpha1)+y_nodes1.*cos(alpha1);
nodes_vect1=[x_nodes1,y_nodes1];

%% Creazione profilo 2
i=2;
[x_nodes2,y_nodes2]=createProfile(CodiceProfilo{i},N_pann2,Chord(i));
x_nodes2=x_nodes2.*cos(alpha2)-y_nodes2.*sin(alpha2);
y_nodes2=x_nodes2.*sin(alpha2)+y_nodes2.*cos(alpha2);
x_nodes2 = x_nodes2 + x12; 
y_nodes2= y_nodes2 + y12;
nodes_vect2=[x_nodes2,y_nodes2];

%% Trovo grandezze importanti dei profili
% Control points profile 1
x_center1 = 0.5*(x_nodes1(1:N_pann1) + x_nodes1(2:N_pann1+1)); 
y_center1 = 0.5*(y_nodes1(1:N_pann1) + y_nodes1(2:N_pann1+1));
centre_vect1=[x_center1 y_center1];
% Control points profile 2
x_center2 = 0.5*(x_nodes2(1:N_pann2) + x_nodes2(2:N_pann2+1)); 
y_center2 = 0.5*(y_nodes2(1:N_pann2) + y_nodes2(2:N_pann2+1));
centre_vect2=[x_center2 y_center2];
% Profile 1 vectors
pan_vect1=zeros(N_pann1,2);
tangent_versor1=zeros(N_pann1,2);
normal_versor1=zeros(N_pann1,2);
lungh_vect1=zeros(N_pann1,1);
for i=1:N_pann1
    pan_vect1(i,:)=nodes_vect1(i+1,:)-nodes_vect1(i,:);
    tangent_versor1(i,:)=pan_vect1(i,:)/norm(pan_vect1(i,:));
    normal_versor1(i,1)=-tangent_versor1(i,2);
    normal_versor1(i,2)=tangent_versor1(i,1);
    lungh_vect1(i,1)=norm(pan_vect1(i,:));
end
% Profile 2 vectors
pan_vect2=zeros(N_pann2,2);
tangent_versor2=zeros(N_pann2,2);
normal_versor2=zeros(N_pann2,2);
lungh_vect2=zeros(N_pann2,1);
for i=1:N_pann2
    pan_vect2(i,:)=nodes_vect2(i+1,:)-nodes_vect2(i,:);
    tangent_versor2(i,:)=pan_vect2(i,:)/norm(pan_vect2(i,:));
    normal_versor2(i,1)=-tangent_versor2(i,2);
    normal_versor2(i,2)=tangent_versor2(i,1);
    lungh_vect2(i,1)=norm(pan_vect2(i,:));
end


figure(); hold on; box on;
plot(x_nodes1,y_nodes1,'k')
plot(x_nodes1(1:10:end),y_nodes1(1:10:end),'ks','MarkerFaceColor','r')
plot(x_nodes2,y_nodes2,'b')
plot(x_nodes2(1:10:end),y_nodes2(1:10:end),'ks','MarkerFaceColor','b')
xlim([0 2.5])
ylim([-1 1])
plot([0 2.5],[0 0],'k','LineWidth',2)
axis off

axis equal

%% Matrice di rotazione profilo 1 (sist locale--> sist globale)
L2G_1=zeros(2*N_pann1,2);
for j=1:N_pann1
    L2G_1((2*j-1:2*j),1)=tangent_versor1(j,:)';
    L2G_1((2*j-1:2*j),2)=normal_versor1(j,:)';
end
% La matrice del pannello i-esimo sarà: R_i=Loc_to_Glob_matrix((i,i+1),:)

%% Matrice di rotazione profilo 2 (sist locale--> sist globale)
L2G_2=zeros(2*N_pann2,2);
for j=1:N_pann2
    L2G_2((2*j-1:2*j),1)=tangent_versor2(j,:)';
    L2G_2((2*j-1:2*j),2)=normal_versor2(j,:)';
end

%% Build the linear system Ax=b
A11 = zeros(N_pann1,N_pann1); 
A22 = zeros(N_pann2,N_pann2);
A12 = zeros(N_pann1,N_pann2); 
A21 = zeros(N_pann2,N_pann1);
a11 = zeros(N_pann1,1); 
a22 = zeros(N_pann2,1);
a12 = zeros(N_pann1,1); 
a21 = zeros(N_pann2,1);
C11 = zeros(1,N_pann1); 
C22 = zeros(1,N_pann2);
C12 = zeros(1,N_pann2); 
C21 = zeros(1,N_pann1);
c11 = 0; 
c22=0; 
c12 = 0; 
c21=0;
b=zeros(N_pann1+N_pann2+2,1);

%% Trovo matrice quadrata A11 e ultima riga di autoinduzione a11
% Autoinduzione profilo 1
for i = 1: N_pann1
    for j = 1: N_pann1
        R_j=L2G_1((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect1(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect1(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        A11(i,j) = dot(Us_j_ri,normal_versor1(i,:));
        a11(i,1) = a11(i,1) + dot(Uv_j_ri,normal_versor1(i,:));
    end
end

%% Trovo matrice quadrata A22 e ultima riga di autoinduzione a22
% Autoinduzione profilo 2
for i = 1: N_pann2
    for j = 1: N_pann2
        R_j=L2G_2((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect2(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect2(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        A22(i,j) = dot(Us_j_ri,normal_versor2(i,:));
        a22(i,1) = a22(i,1) + dot(Uv_j_ri,normal_versor2(i,:));
    end
end

%% Trovo A12 e a12
% Induzione profilo 2 su profilo 1
for i = 1: N_pann1
    for j = 1: N_pann2
        R_j=L2G_2((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect1(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect1(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        A12(i,j) = dot(Us_j_ri,normal_versor1(i,:));
        a12(i,1) = a12(i,1) + dot(Uv_j_ri,normal_versor1(i,:));
    end
end

%% Trovo A21 e a21
% Induzione profilo 1 su profilo 2
for i = 1: N_pann2
    for j = 1: N_pann1
        R_j=L2G_1((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect2(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect2(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        A21(i,j) = dot(Us_j_ri,normal_versor2(i,:));
        a21(i,1) = a21(i,1) + dot(Uv_j_ri,normal_versor2(i,:));
    end
end

%% Costruisco C11 and c11 (condizione di Kutta)
% Autoinduzione profilo 1
for j = 1:N_pann1
    R_j=L2G_1((2*j-1:2*j),:);
    [Us_j_ri] = ViSorgente(centre_vect1(1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect1(1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    C11(1,j)=dot(Us_j_ri,tangent_versor1(1,:));
    c11=c11+dot(Uv_j_ri,tangent_versor1(1,:));
    [Us_j_ri] = ViSorgente(centre_vect1(N_pann1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect1(N_pann1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    C11(1,j)=C11(1,j)+dot(Us_j_ri,tangent_versor1(N_pann1,:));
    c11=c11+dot(Uv_j_ri,tangent_versor1(N_pann1,:));   
end

%% Costruisco C22 and c22 (condizione di Kutta)
% Autoinduzione profilo 2
for j = 1:N_pann2
    R_j=L2G_2((2*j-1:2*j),:);
    [Us_j_ri] = ViSorgente(centre_vect2(1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect2(1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    C22(1,j)=dot(Us_j_ri,tangent_versor2(1,:));
    c22=c22+dot(Uv_j_ri,tangent_versor2(1,:));
    [Us_j_ri] = ViSorgente(centre_vect2(N_pann2,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect2(N_pann2,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    C22(1,j)=C22(1,j)+dot(Us_j_ri,tangent_versor2(N_pann2,:));
    c22=c22+dot(Uv_j_ri,tangent_versor2(N_pann2,:));    
end

%% build C12 and c12
% Induzione profilo 2 su profilo 1
for j = 1:N_pann2
    R_j=L2G_2((2*j-1:2*j),:);
    [Us_j_ri] = ViSorgente(centre_vect1(1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect1(1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    C12(1,j)=dot(Us_j_ri,tangent_versor1(1,:));
    c12=c12+dot(Uv_j_ri,tangent_versor1(1,:));
    [Us_j_ri] = ViSorgente(centre_vect1(N_pann1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect1(N_pann1,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
    C12(1,j)=C12(1,j)+dot(Us_j_ri,tangent_versor1(N_pann1,:));
    c12=c12+dot(Uv_j_ri,tangent_versor1(N_pann1,:));   
end

%% build C21 and c21
% Induzione profilo 1 su profilo 2
for j = 1:N_pann1
    R_j=L2G_1((2*j-1:2*j),:);
    [Us_j_ri] = ViSorgente(centre_vect2(1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect2(1,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    C21(1,j)=dot(Us_j_ri,tangent_versor2(1,:));
    c21=c21+dot(Uv_j_ri,tangent_versor2(1,:));
    [Us_j_ri] = ViSorgente(centre_vect2(N_pann2,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    [Uv_j_ri] = ViVortice(centre_vect2(N_pann2,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
    C21(1,j)=C21(1,j)+dot(Us_j_ri,tangent_versor2(N_pann2,:));
    c21=c21+dot(Uv_j_ri,tangent_versor2(N_pann2,:));  
end

%% Assemblo la matrice A
A=[A11 A12 a11 a12;
   A21 A22 a21 a22;
   C11 C12 c11 c12;
   C21 C22 c21 c22];

%% Scrivo termine noto sistema lineare
for i = 1: N_pann1
    b(i) = -dot(U_inf,normal_versor1(i,:));
end
for i = 1: N_pann2
    b(N_pann1+i) = -dot(U_inf,normal_versor2(i,:));
end
b(N_pann1+N_pann2+1)=-dot(U_inf,tangent_versor1(1,:))-dot(U_inf,tangent_versor1(N_pann1,:));
b(N_pann1+N_pann2+2)=-dot(U_inf,tangent_versor2(1,:))-dot(U_inf,tangent_versor2(N_pann2,:));

%% Risolvo il sistema lineare
x_vect=linsolve(A,b);

sigma1 = x_vect(1:N_pann1);
sigma2  =x_vect(N_pann1+1:N_pann1+N_pann2);
gamma1= x_vect(N_pann1+N_pann2+1);
gamma2 = x_vect(N_pann1+N_pann2+2);

%% Calcolo velocità nei punti di controllo
u_ri_1= zeros(N_pann1,1); 
v_ri_1 = u_ri_1;
u_ri_2= zeros(N_pann2,1); 
v_ri_2 = u_ri_2;
% airfoil 1
for i =1:N_pann1
    u_ri_1(i) = U_inf(1); 
    v_ri_1(i)=U_inf(2);
    for j = 1:N_pann1
        R_j=L2G_1((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect1(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect1(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        u_ri_1(i) = u_ri_1(i) + sigma1(j)*Us_j_ri(1,1) + gamma1*Uv_j_ri(1,1);
        v_ri_1(i) = v_ri_1(i) + sigma1(j)*Us_j_ri(2,1) + gamma1*Uv_j_ri(2,1);        
    end
    for j = 1:N_pann2
        R_j=L2G_2((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect1(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect1(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
       u_ri_1(i) = u_ri_1(i) + sigma2(j)*Us_j_ri(1,1) + gamma2*Uv_j_ri(1,1);
       v_ri_1(i) = v_ri_1(i) + sigma2(j)*Us_j_ri(2,1) + gamma2*Uv_j_ri(2,1);        
    end
end
% arifoil 2
for i =1:N_pann2
    u_ri_2(i) = U_inf(1); 
    v_ri_2(i)=U_inf(2);
    for j = 1:N_pann1
        R_j=L2G_1((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect2(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect2(i,:)',nodes_vect1(j,:)',...
            nodes_vect1(j+1,:)',R_j,R_j');
       u_ri_2(i) = u_ri_2(i) + sigma1(j)* Us_j_ri(1,1)+ gamma1*Uv_j_ri(1,1);
       v_ri_2(i) = v_ri_2(i) + sigma1(j)* Us_j_ri(2,1)+ gamma1*Uv_j_ri(2,1);
    end
    for j = 1:N_pann2
        R_j=L2G_2((2*j-1:2*j),:);
        [Us_j_ri] = ViSorgente(centre_vect2(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
        [Uv_j_ri] = ViVortice(centre_vect2(i,:)',nodes_vect2(j,:)',...
            nodes_vect2(j+1,:)',R_j,R_j');
       u_ri_2(i) = u_ri_2(i) + sigma2(j)*Us_j_ri(1,1) + gamma2*Uv_j_ri(1,1);
       v_ri_2(i) = v_ri_2(i) + sigma2(j)*Us_j_ri(2,1)+ gamma2*Uv_j_ri(2,1);       
    end
end

if (max(u_ri_1.*normal_versor1(:,1) + v_ri_1.*normal_versor1(:,2))>10^(-14)) ...
        || (max(u_ri_2.*normal_versor2(:,1) + v_ri_2.*normal_versor2(:,2))>10^(-14))
    disp('There is a bug in the program!')
    stop
end
U_ri1=[u_ri_1 v_ri_1];
U_ri2=[u_ri_2 v_ri_2];

%% Calcolo Cp sui 2 profili
Cp1=zeros(N_pann1,1);
Cp2=zeros(N_pann2,1);
for i=1:N_pann1
    Cp1(i,1)=1-norm(U_ri1(i,:))^2/norm(U_inf)^2;
end
for i=1:N_pann2
    Cp2(i,1)=1-norm(U_ri2(i,:))^2/norm(U_inf)^2;
end

%% Calcolo Cl sui 2 profili
Cl1=0;
for i=1:N_pann1
    Chord1=Chord(1);
    Cl1=Cl1-dot(((1/Chord1)*Cp1(i,1)*lungh_vect1(i).*normal_versor1(i,:)),[0;1]);
end
Cl2=0;
for i=1:N_pann2
    Chord2=Chord(2);
    Cl2=Cl2-dot(((1/Chord2)*Cp2(i,1)*lungh_vect2(i).*normal_versor2(i,:)),[0;1]);
end


