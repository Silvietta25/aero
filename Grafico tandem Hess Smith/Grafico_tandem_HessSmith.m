%% Grafico tandem Hess-Smith
clear all
close all
addpath mat_functions
clc

%% Densità
rho=1.225;

%% Input
U = 1;
alpha1 = 5*pi/180;
U_inf(1)=U; 
U_inf(2)=0;

delta_alpha_vect=0:0.5:7;
alpha2_vect = alpha1+delta_alpha_vect.*pi/180;
Cl_tot_vect=[];

for l=1:length(alpha2_vect)
    alpha2=alpha2_vect(l);
    %% Dati profili
    TestCase = 0;
    NCorpi = 2;  % Numero di corpi da analizzare
    CodiceProfilo = cell(NCorpi, 1);
    CodiceProfilo{1} = '0024';
    CodiceProfilo{2} = '0012';
    Chord = [1 1];
    N_pann1 = 101;
    N_pann2 = 101;
    
    %% Creazione profilo 1
    i=1;
    [x_nodes1,y_nodes1]=createProfile(CodiceProfilo{i},N_pann1,Chord(i));
    x_nodes1=(x_nodes1).*cos(alpha1)+y_nodes1.*sin(alpha1);
    y_nodes1=-(x_nodes1).*sin(alpha1)+y_nodes1.*cos(alpha1);
    nodes_vect1=[x_nodes1,y_nodes1];
    
    %% Distanza tra i profili in tandem
    Gap=0.05;
    x12=x_nodes1(1)+Gap;
    y12=y_nodes1(1);
    
    %% Creazione profilo 2
    i=2;
    [x_nodes2,y_nodes2]=createProfile(CodiceProfilo{i},N_pann2,Chord(i));
    x_nodes2=(x_nodes2).*cos(alpha2)+y_nodes2.*sin(alpha2);
    y_nodes2=-(x_nodes2).*sin(alpha2)+y_nodes2.*cos(alpha2);
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
    
    %% Plot del Cp dei profili
    % figure
    % plot(centre_vect1((1:(N_pann1+1)/2),1),-Cp1((1:(N_pann1+1)/2),1),'r','LineWidth',2)
    % grid on
    % xlim([0,1])
    % hold on
    % plot(centre_vect1((N_pann1+1)/2:end,1),-Cp1((N_pann1+1)/2:end,1),'b','LineWidth',2)
    % hold on
    % plot(centre_vect2((1:(N_pann2+1)/2),1)-x12,-Cp2((1:(N_pann2+1)/2),1),'g','LineWidth',2)
    % hold on
    % plot(centre_vect2((N_pann2+1)/2:end,1)-x12,-Cp2((N_pann2+1)/2:end,1),'y','LineWidth',2)
    % legend('-Cp ventre profilo 1','-Cp dorso profilo 1','-Cp ventre profilo 2','-Cp dorso profilo 2')
    % 
    %% Calcolo Cl sui 2 profili: adimensionalizzo con la corda maggiore tra le due
    Gamma1=sum(gamma1.*lungh_vect1);
    L1_2D=rho*norm(U_inf)*Gamma1;
    Cl1_KJ=2*Gamma1/(norm(U_inf)*Chord(1));
    Cl1=0;
    for i=1:N_pann1
        Cl1=Cl1-dot((1/Chord(1))*Cp1(i,1)*lungh_vect1(i).*normal_versor1(i,:),[0;1/cos(alpha1)]);
    end
    Gamma2=sum(gamma2.*lungh_vect2);
    L2_2D=rho*norm(U_inf)*Gamma2;
    Cl2_KJ=2*Gamma2/(norm(U_inf)*Chord(2));
    Cl2=0;
    for i=1:N_pann2
        Cl2=Cl2-dot(((1/Chord(2))*Cp2(i,1)*lungh_vect2(i).*normal_versor2(i,:)),[0;1/cos(alpha2)]);
    end
    
    Cl_tot_vect(end+1)=Cl1;
end

%% Plotto il Cl al variare di alpha2 (angolo assoluto)
% Come Cl del tandem si fa riferimento al Cl del corpo 1
% L'angolo dell'ala è fissato ad incidenza di 8°
figure
plot(alpha2_vect*180/pi,Cl_tot_vect,'LineWidth',5,'Color','b')
hold on
yline(0.6602,'LineWidth',5,'Color','g','LineStyle','--')
grid on
legend('Cl tandem','Cl senza tandem','FontSize',40)
xlabel('\alpha_2','FontSize',40)
ylabel('C_L','FontSize',40)
xlim([5 12])
ylim([0 3.5])
ax=gca;
ax.FontSize=40;
title('Cl al variare di \alpha_2', 'FontSize',50)


