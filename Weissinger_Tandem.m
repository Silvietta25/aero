close all
clear all
clc

%% Test Case 1

U_Inf_Mag = 1;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 2; % numero corpi

config.RootChord = [1 0.5]; % corde delle 2 ali
config.DihedralAngle = [0 0]; % angolo di diedro in gradi delle 2 ali
config.SweepAngle = [10 10]; % angolo di freccia in gradi delle 2 ali
config.TaperRatio = [1 0.5]; % rapporto di rastremazione delle 2 ali
config.Span = [2 1]; % apertura alare delle 2 ali
config.LEPosition_X = [0 5]; % posizione su x dell'origine delle 2 ali
config.LEPosition_Y = [0 0]; % posizione su y dell'origine delle 2 ali
config.LEPosition_Z = [0 1]; % posizione su z dell'origine delle 2 ali

config.RotationAngle_X = [0 0]; % angolo di rollio (lo lascio sempre nullo)
config.RotationAngle_Y = [10 25]; % incidenza delle 2 ali
config.RotationAngle_Z = [0 0]; % angolo d'imbardata (lo lascio sempre nullo)

% Discretization options
config.SemiSpanwiseDiscr = [5 5]; % discretizzazione su semiapertura alare per le 2 ali
config.ChordwiseDiscr = [5 5]; % discretizzazione su corda per le 2 ali


%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);


for iCorpo = 1:config.NCorpi

    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);

end

%% Trovo punti ad 1/4 della corda del Corpo 1
QuarterPoints_dx_Corpo1=cell(1,config.SemiSpanwiseDiscr(1));
TipQuarter_Corpo1=WingExtremes{1,1}{3,1}.LE-(WingExtremes{1,1}{3,1}.LE-WingExtremes{1,1}{3,1}.TE)./4;
RootQuarter_Corpo1=WingExtremes{1,1}{2,1}.LE-(WingExtremes{1,1}{2,1}.LE-WingExtremes{1,1}{2,1}.TE)./4;
TipQuarter_sx_Corpo1=[TipQuarter_Corpo1(1) -TipQuarter_Corpo1(2) TipQuarter_Corpo1(3)];
LEDiscr_Corpo1 = linspace(1, 0, config.SemiSpanwiseDiscr(1)+1);
AerodinCentre_Corpo1 = @(x) (TipQuarter_Corpo1-RootQuarter_Corpo1).*x+RootQuarter_Corpo1;
for j=2:length(LEDiscr_Corpo1)
    QuarterPoints_dx_Corpo1{1,j-1}=AerodinCentre_Corpo1((LEDiscr_Corpo1(j)+LEDiscr_Corpo1(j-1))/2);
end
QuarterPoints_sx_Corpo1=flip(QuarterPoints_dx_Corpo1);
for j=1:length(QuarterPoints_sx_Corpo1)
    QuarterPoints_sx_Corpo1{1,j}(2)=-QuarterPoints_sx_Corpo1{1,j}(2);
end
QuarterPoints_Corpo1=[QuarterPoints_dx_Corpo1,QuarterPoints_sx_Corpo1];

%% Trovo punti ad 1/4 della corda del Corpo 2
QuarterPoints_dx_Corpo2=cell(1,config.SemiSpanwiseDiscr(2));
TipQuarter_Corpo2=WingExtremes{2,1}{3,1}.LE-(WingExtremes{2,1}{3,1}.LE-WingExtremes{2,1}{3,1}.TE)./4;
RootQuarter_Corpo2=WingExtremes{2,1}{2,1}.LE-(WingExtremes{2,1}{2,1}.LE-WingExtremes{2,1}{2,1}.TE)./4;
TipQuarter_sx_Corpo2=[TipQuarter_Corpo2(1) -TipQuarter_Corpo2(2) TipQuarter_Corpo2(3)];
LEDiscr_Corpo2 = linspace(1, 0, config.SemiSpanwiseDiscr(2)+1);
AerodinCentre_Corpo2 = @(x) (TipQuarter_Corpo2-RootQuarter_Corpo2).*x+RootQuarter_Corpo2;
for j=2:length(LEDiscr_Corpo2)
    QuarterPoints_dx_Corpo2{1,j-1}=AerodinCentre_Corpo2((LEDiscr_Corpo1(j)+LEDiscr_Corpo1(j-1))/2);
end
QuarterPoints_sx_Corpo2=flip(QuarterPoints_dx_Corpo2);
for j=1:length(QuarterPoints_sx_Corpo2)
    QuarterPoints_sx_Corpo2{1,j}(2)=-QuarterPoints_sx_Corpo2{1,j}(2);
end
QuarterPoints_Corpo2=[QuarterPoints_dx_Corpo2,QuarterPoints_sx_Corpo2];
QuarterPoints{1,1}=QuarterPoints_Corpo1;
QuarterPoints{2,1}=QuarterPoints_Corpo2;

%% Plot Corpi e centri aerodinamici
figure
for numCorpo=1:config.NCorpi
    for i=1:config.ChordwiseDiscr(numCorpo)
        for j=1:2*config.SemiSpanwiseDiscr(numCorpo)
            x_vect=[internalMesh{numCorpo,1}{i,j}.LERoot(1); internalMesh{numCorpo,1}{i,j}.LEtip(1); ...
                    internalMesh{numCorpo,1}{i,j}.TEtip(1); internalMesh{numCorpo,1}{i,j}.TERoot(1)];
            y_vect=[internalMesh{numCorpo,1}{i,j}.LERoot(2); internalMesh{numCorpo,1}{i,j}.LEtip(2); ...
                    internalMesh{numCorpo,1}{i,j}.TEtip(2); internalMesh{numCorpo,1}{i,j}.TERoot(2)];
            z_vect=[internalMesh{numCorpo,1}{i,j}.LERoot(3); internalMesh{numCorpo,1}{i,j}.LEtip(3); ...
                    internalMesh{numCorpo,1}{i,j}.TEtip(3); internalMesh{numCorpo,1}{i,j}.TERoot(3)];
            plot3(x_vect,y_vect,z_vect,'r')
            axis equal
            hold on
            grid on
            x_c=ControlPoints{numCorpo,1}{i,j}.Coords(1);
            y_c=ControlPoints{numCorpo,1}{i,j}.Coords(2);
            z_c=ControlPoints{numCorpo,1}{i,j}.Coords(3);
            plot3(x_c,y_c,z_c,'.','MarkerSize',8,'Color','g')
            hold on
        end
    end
    hold on
    for h=1:2*config.SemiSpanwiseDiscr(numCorpo)
        x_q=QuarterPoints{numCorpo,1}{1,h}(1);
        y_q=QuarterPoints{numCorpo,1}{1,h}(2);
        z_q=QuarterPoints{numCorpo,1}{1,h}(3);
        plot3(x_q,y_q,z_q,'.','MarkerSize',8,'Color','b')
        hold on
    end
end

%% Matrices initialization
NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            QuarterPoint = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            
            for jCorpo = 1:config.NCorpi
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(QuarterPoint, Extreme_1, Extreme_2);

                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(QuarterPoint, Extreme_1, Extreme_2);
                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(QuarterPoint, Extreme_1, Extreme_2);
                        
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                      
                    end
                end
            end   
        end
    end
end

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            % Update row index
            rowIndex = rowIndex + 1;
  
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
        end
    end
end

%% Solve the linear system
Solution = linsolve(matriceA, TermineNoto);
Gamma = cell(config.NCorpi, 1);

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end 
    end 
end
%% Calcolo Lift 2D e 3D Corpo 1 (quello anteriore)
Gamma_i_Corpo1=[];
for i=1:2*config.SemiSpanwiseDiscr(1)
    Gamma_i_Corpo1(i)=sum(Gamma{1,1}(:,i));
end
L_2D_i_Corpo1=rho*norm(U_Inf).*Gamma_i_Corpo1';
N1=config.ChordwiseDiscr(1);
M1=config.SemiSpanwiseDiscr(1);
% Calcolo il Cl complessivo (adimensionalizzo con apertura alare e media
% algebrica corde del Corpo 1)
delta_b_vect1=zeros(2*M1,1);
for i=1:2*M1
    delta_b_vect1(i)=abs(internalMesh{1,1}{1,i}.LEtip(2)-internalMesh{1,1}{1,i}.LERoot(2));
end
L_3D_i_Corpo1=L_2D_i_Corpo1.*delta_b_vect1;
L_3D_Corpo1=sum(L_3D_i_Corpo1);
c_media1=(config.TipChord(1)+config.RootChord(1))/2;
Cl_3D_Corpo1=L_3D_Corpo1/(0.5*rho*(norm(U_Inf)^2)*config.Span(1)*c_media1);

%% Calcolo Lift 2D e 3D Corpo 2 (quello anteriore)
Gamma_i_Corpo2=[];
for i=1:2*config.SemiSpanwiseDiscr(2)
    Gamma_i_Corpo2(i)=sum(Gamma{2,1}(:,i));
end
L_2D_i_Corpo2=rho*norm(U_Inf).*Gamma_i_Corpo2';
N2=config.ChordwiseDiscr(2);
M2=config.SemiSpanwiseDiscr(2);
% Calcolo il Cl complessivo (adimensionalizzo con apertura alare e media
% algebrica corde del Corpo 2)
delta_b_vect2=zeros(2*M2,1);
for i=1:2*M2
    delta_b_vect2(i)=abs(internalMesh{2,1}{1,i}.LEtip(2)-internalMesh{2,1}{1,i}.LERoot(2));
end
L_3D_i_Corpo2=L_2D_i_Corpo2.*delta_b_vect2;
L_3D_Corpo2=sum(L_3D_i_Corpo2);
c_media2=(config.TipChord(2)+config.RootChord(2))/2;
Cl_3D_Corpo2=L_3D_Corpo2/(0.5*rho*(norm(U_Inf)^2)*config.Span(2)*c_media2);

%% Calcolo velocità indotta da scia su punto ad 1/4 della sezione i-esima
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
        % Update row index
        rowIndex = rowIndex + 1;
        columnIndex = 0;
        
        QuarterPoint = QuarterPoints{iCorpo,1}{1,SpanPanel_i};
        
        for jCorpo = 1:config.NCorpi
            
            % Cycle on all of its chordwise panels
            for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                % Cycle on all of its spanwise panels
                for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                    
                    % Update column index
                    columnIndex = columnIndex + 1;
                    
                    % Compute the influence induced by first
                    % semi-infinite vortex
                    Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                    Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                    U = vortexInfluence(QuarterPoint, Extreme_1, Extreme_2);

                    % Compute the influence induced by second
                    % semi-infinite vortex
                    Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                    Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                    U = U + vortexInfluence(QuarterPoint, Extreme_1, Extreme_2);
                    q_jk_wake{columnIndex,rowIndex}=U;

                end
            end
        end  
    end
end

%% Calcolo velocità indotta da scia su sezione i-esima del corpo 1
v_ind_i_Corpo1=zeros(2*M1,3); % velocità indotta da scia su centro aerodinamico 
% (cioè ad 1/4 della corda) del profilo i-esimo (riga=posizione su
% apertura, colonna=componente della velocità)
for j=1:NPanelsTot
    for k=1:2*config.SemiSpanwiseDiscr(1)
    v_ind_jk_wake_Corpo1{j,k}=q_jk_wake{j,k}.*Solution(j,1); % velocità indotta 
    % dal pannello generico (j-esimo) sul centro aerodinamico della sezione k-esima
    end
end
% Trovo velocità indotta sulla sezione i-esima: v_ind_i
for i=1:2*M1
    for j=1:NPanelsTot
        v_ind_i_Corpo1(i,:)=v_ind_i_Corpo1(i,:)+v_ind_jk_wake_Corpo1{j,i};
    end
end

%% Calcolo velocità indotta da scia su sezione i-esima del corpo 1
v_ind_i_Corpo2=zeros(2*M2,3); % velocità indotta da scia su centro aerodinamico 
% (cioè ad 1/4 della corda) del profilo i-esimo (riga=posizione su
% apertura, colonna=componente della velocità)
for j=1:NPanelsTot
    for k=1:2*config.SemiSpanwiseDiscr(2)
    v_ind_jk_wake_Corpo2{j,k}=q_jk_wake{j,2*config.SemiSpanwiseDiscr(2)+k}.*Solution(j,1); % velocità indotta 
    % dal pannello generico (j-esimo) sul centro aerodinamico della sezione k-esima
    end
end
% Trovo velocità indotta sulla sezione i-esima: v_ind_i
for i=1:2*M2
    for j=1:NPanelsTot
        v_ind_i_Corpo2(i,:)=v_ind_i_Corpo2(i,:)+v_ind_jk_wake_Corpo2{j,i};
    end
end

%% Calcolo incidenza indotta da scia su sezione i-esima del Corpo 1
alpha_ind_i_Corpo1=zeros(2*M1,1);
for i=1:2*M1 % ciclo su apertura alare (da dx a sx)
    alpha_ind_i_Corpo1(i,1)=atand(v_ind_i_Corpo1(i,3)*Normals{1,1}{1,1}.Coords(3)/norm(U_Inf));
end

%% Calcolo incidenza indotta da scia su sezione i-esima del Corpo 2
alpha_ind_i_Corpo2=zeros(2*M2,1);
for i=1:2*M2 % ciclo su apertura alare (da dx a sx)
    alpha_ind_i_Corpo2(i,1)=atand(v_ind_i_Corpo2(i,3)*Normals{2,1}{1,1}.Coords(3)/norm(U_Inf));
end

%% Calcolo resistenza indotta 2D e 3D del Corpo 1
D_2D_Corpo1=zeros(2*M1,1);
for i=1:2*M1
    D_2D_Corpo1(i,1)=L_2D_i_Corpo1(i)*sind(abs(alpha_ind_i_Corpo1(i)));
end
Cd_2D_Corpo1=D_2D_Corpo1./(0.5*rho*(norm(U_Inf)^2)*c_media1);
D_3D_i_Corpo1=D_2D_Corpo1.*delta_b_vect1;
D_3D_Corpo1=sum(D_3D_i_Corpo1);
Cd_3D_Corpo1=D_3D_Corpo1/(0.5*rho*(norm(U_Inf)^2)*config.Span(1)*c_media1);

%% Calcolo resistenza indotta 2D e 3D del Corpo 2
D_2D_Corpo2=zeros(2*M2,1);
for i=1:2*M2
    D_2D_Corpo2(i,1)=L_2D_i_Corpo2(i)*sind(abs(alpha_ind_i_Corpo2(i)));
end
Cd_2D_Corpo2=D_2D_Corpo2./(0.5*rho*(norm(U_Inf)^2)*c_media2);
D_3D_i_Corpo2=D_2D_Corpo2.*delta_b_vect2;
D_3D_Corpo2=sum(D_3D_i_Corpo2);
Cd_3D_Corpo2=D_3D_Corpo2/(0.5*rho*(norm(U_Inf)^2)*config.Span(2)*c_media2);



