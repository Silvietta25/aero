close all
clear all
clc

%% Test Case 1

U_Inf_Mag = 1;
beta = 0;
U_Inf = [cosd(beta) -sind(beta) 0] .* U_Inf_Mag; % cambio segno al seno perchè
% voglio beta positivo se U_inf arriva da lato della semiala destra
rho = 1.225;

config.NCorpi = 1;

config.RootChord = [10];
config.DihedralAngle = [5]; % [°]
config.SweepAngle = [10]; % [°]
config.TaperRatio = [0.6];
config.Span = 80;
config.LEPosition_X = [0];
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];

config.RotationAngle_X = [0]; % per dare rollio (ma devo modificare 
% createStructure per usarlo)
ALPHA=5;
config.RotationAngle_Y = [ALPHA]; % per dare calettamento
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [20]; % parti in cui suddivido la semi-apertura alare (M pannelli)
config.ChordwiseDiscr = [10]; % parti in cui suddivido la corda (N pannelli)


%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Corda media aerodinamica dell'ala (utile per adimensionalizzare poi)
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
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, ...
        InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, ...
        WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end

%% Trovo centro aerodinamico delle varie sezioni
QuarterPoints_dx=cell(1,config.SemiSpanwiseDiscr);
TipQuarter=WingExtremes{1,1}{3,1}.LE-(WingExtremes{1,1}{3,1}.LE-WingExtremes{1,1}{3,1}.TE)./4;
RootQuarter=WingExtremes{1,1}{2,1}.LE-(WingExtremes{1,1}{2,1}.LE-WingExtremes{1,1}{2,1}.TE)./4;
TipQuarter_sx=[TipQuarter(1) -TipQuarter(2) TipQuarter(3)];
LEDiscr = linspace(1, 0, config.SemiSpanwiseDiscr(iCorpo)+1);
AerodinCentre = @(x) (TipQuarter-RootQuarter).*x+RootQuarter;
for j=2:length(LEDiscr)
    QuarterPoints_dx{1,j-1}=AerodinCentre((LEDiscr(j)+LEDiscr(j-1))/2);
end
QuarterPoints_sx=flip(QuarterPoints_dx);
for j=1:length(QuarterPoints_sx)
    QuarterPoints_sx{1,j}(2)=-QuarterPoints_sx{1,j}(2);
end
QuarterPoints=[QuarterPoints_dx,QuarterPoints_sx];

%% Plot pianta alare
figure
for i=1:config.ChordwiseDiscr
    for j=1:2*config.SemiSpanwiseDiscr
        x_vect=[internalMesh{1,1}{i,j}.LERoot(1); internalMesh{1,1}{i,j}.LEtip(1); ...
                internalMesh{1,1}{i,j}.TEtip(1); internalMesh{1,1}{i,j}.TERoot(1)];
        y_vect=[internalMesh{1,1}{i,j}.LERoot(2); internalMesh{1,1}{i,j}.LEtip(2); ...
                internalMesh{1,1}{i,j}.TEtip(2); internalMesh{1,1}{i,j}.TERoot(2)];
        z_vect=[internalMesh{1,1}{i,j}.LERoot(3); internalMesh{1,1}{i,j}.LEtip(3); ...
                internalMesh{1,1}{i,j}.TEtip(3); internalMesh{1,1}{i,j}.TERoot(3)];
        plot3(x_vect,y_vect,z_vect,'r')
        axis equal
        hold on
        grid on
        x_c=ControlPoints{1,1}{i,j}.Coords(1);
        y_c=ControlPoints{1,1}{i,j}.Coords(2);
        z_c=ControlPoints{1,1}{i,j}.Coords(3);
        plot3(x_c,y_c,z_c,'.','MarkerSize',8,'Color','g')
        hold on
    end
end
hold on
for j=1:2*config.SemiSpanwiseDiscr
    x_q=QuarterPoints{1,j}(1);
    y_q=QuarterPoints{1,j}(2);
    z_q=QuarterPoints{1,j}(3);
    plot3(x_q,y_q,z_q,'.','MarkerSize',8,'Color','b')
    hold on
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
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
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
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);


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
Gamma_i=[];
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

%% Compute the 2D and 3D Lift
for i=1:2*config.SemiSpanwiseDiscr
    Gamma_i(i)=sum(Gamma{1,1}(:,i));
end
L_2D=rho*norm(U_Inf).*Gamma_i';
N=config.ChordwiseDiscr; % numero di pannelli rispetto a discretizzazione 
% in direzione della corda
M=config.SemiSpanwiseDiscr; % numero di pannelli rispetto a discretizzazione
% su semiapertura alare
% Adimensionalizzo e trovo il Cl_2D
Cl_2D=L_2D./(0.5*rho*(norm(U_Inf)^2)*config.MAC);
% Calcolo il Cl complessivo
delta_b_vect=zeros(2*M,1);
for i=1:2*M
    delta_b_vect(i)=abs(internalMesh{1,1}{1,i}.LEtip(2)-internalMesh{1,1}{1,i}.LERoot(2));
end
L_3D_i=L_2D.*delta_b_vect;
L_3D=sum(L_3D_i);
c_media=(config.TipChord+config.RootChord)/2;
Cl_3D=L_3D/(0.5*rho*(norm(U_Inf)^2)*config.Span*c_media);

%% Calcolo velocità indotta da scia su punto ad 1/4 della sezione i-esima
rowIndex = 0;
% q_jk_wake=cell(NPanelsTot,2*config.SemiSpanwiseDiscr(iCorpo));
for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
    % Update row index
    rowIndex = rowIndex + 1;
    columnIndex = 0;
    
    QuarterPointHere = QuarterPoints{1,SpanPanel_i};
        
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
            U = vortexInfluence(QuarterPointHere, Extreme_1, Extreme_2);

            % Compute the influence induced by second
            % semi-infinite vortex
            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
            U = U + vortexInfluence(QuarterPointHere, Extreme_1, Extreme_2);
            q_jk_wake{columnIndex,rowIndex}=U;
        end
    end        
end

%% Calcolo velocità indotta da scia su sezione i-esima
v_ind_i=zeros(2*M,3); % velocità indotta da scia su centro aerodinamico 
% (cioè ad 1/4 della corda) del profilo i-esimo (riga=posizione su
% apertura, colonna=componente della velocità)
for j=1:NPanelsTot
    for k=1:2*config.SemiSpanwiseDiscr(iCorpo)
    v_ind_jk_wake{j,k}=q_jk_wake{j,k}.*Solution(j,1); % velocità indotta 
    % dal pannello generico (j-esimo) sul centro aerodinamico della sezione k-esima
    end
end
% Trovo velocità indotta sulla sezione i-esima: v_ind_i
for i=1:2*M
    for j=1:NPanelsTot
        v_ind_i(i,:)=v_ind_i(i,:)+v_ind_jk_wake{j,i};
    end
end
%% Calcolo incidenza indotta da scia su sezione i-esima
alpha_ind=zeros(2*M,1);
for i=1:2*M % ciclo su apertura alare (da dx a sx)
    alpha_ind(i,1)=atand(v_ind_i(i,3)*Normals{1,1}{1,1}.Coords(3)/norm(U_Inf));
end
%% Calcolo resistenza indotta 2D e 3D
D_2D=zeros(2*M,1);
for i=1:2*M
    D_2D(i,1)=L_2D(i)*sind(abs(alpha_ind(i)));
end
Cd_2D=D_2D./(0.5*rho*(norm(U_Inf)^2)*c_media);
D_3D_i=D_2D.*delta_b_vect;
D_3D=sum(D_3D_i);
Cd_3D=D_3D/(0.5*rho*(norm(U_Inf)^2)*config.Span*c_media);

