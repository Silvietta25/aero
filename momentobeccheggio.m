function [M_CG_tot,Cm_CG_tot]=momentobeccheggio(config,ALPHA1_trim,ALPHA2_trim,Perturbaz_Alpha_velivolo,U_Inf,rho)

ALPHA1=Perturbaz_Alpha_velivolo+ALPHA1_trim;
ALPHA2=Perturbaz_Alpha_velivolo+ALPHA2_trim;
config.RotationAngle_Y = [ALPHA1 ALPHA2];
% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

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


NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

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

alpha_ind_i_Corpo1=zeros(2*M1,1);
for i=1:2*M1 % ciclo su apertura alare (da dx a sx)
alpha_ind_i_Corpo1(i,1)=atand(v_ind_i_Corpo1(i,3)*Normals{1,1}{1,1}.Coords(3)/norm(U_Inf));
end

alpha_ind_i_Corpo2=zeros(2*M2,1);
for i=1:2*M2 % ciclo su apertura alare (da dx a sx)
alpha_ind_i_Corpo2(i,1)=atand(v_ind_i_Corpo2(i,3)*Normals{2,1}{1,1}.Coords(3)/norm(U_Inf));
end

D_2D_Corpo1=zeros(2*M1,1);
for i=1:2*M1
D_2D_Corpo1(i,1)=abs(L_2D_i_Corpo1(i))*sind(abs(alpha_ind_i_Corpo1(i)));
end
Cd_2D_Corpo1=D_2D_Corpo1./(0.5*rho*(norm(U_Inf)^2)*c_media1);
D_3D_i_Corpo1=D_2D_Corpo1.*delta_b_vect1;
D_3D_Corpo1=sum(D_3D_i_Corpo1);
Cd_3D_Corpo1=D_3D_Corpo1/(0.5*rho*(norm(U_Inf)^2)*config.Span(1)*c_media1);

D_2D_Corpo2=zeros(2*M2,1);
for i=1:2*M2
D_2D_Corpo2(i,1)=abs(L_2D_i_Corpo2(i))*sind(abs(alpha_ind_i_Corpo2(i)));
end
Cd_2D_Corpo2=D_2D_Corpo2./(0.5*rho*(norm(U_Inf)^2)*c_media2);
D_3D_i_Corpo2=D_2D_Corpo2.*delta_b_vect2;
D_3D_Corpo2=sum(D_3D_i_Corpo2);
Cd_3D_Corpo2=D_3D_Corpo2/(0.5*rho*(norm(U_Inf)^2)*config.Span(2)*c_media2);

L_tot=L_3D_Corpo1+L_3D_Corpo2;
% Adimensionalizziamo rispetto alle grandezze dell'ala
CL_tot=L_tot/(0.5*rho*(norm(U_Inf)^2*c_media1*2*config.SemiSpan(1)));

D_tot=D_3D_Corpo1+D_3D_Corpo2;
% Adimensionalizziamo rispetto alle grandezze dell'ala
CD_tot=D_tot/(0.5*rho*(norm(U_Inf)^2*c_media1*2*config.SemiSpan(1)));

CA_ala=WingExtremes{1,1}{2,1}.LE+(WingExtremes{1,1}{2,1}.TE-WingExtremes{1,1}{2,1}.LE)/4;
M_CG_ala=-L_3D_Corpo1*(CA_ala(1));
CA_coda=WingExtremes{2,1}{2,1}.LE+(WingExtremes{2,1}{2,1}.TE-WingExtremes{2,1}{2,1}.LE)/4;
M_CG_coda=-L_3D_Corpo2*((CA_coda(1)));
M_CG_tot=M_CG_ala+M_CG_coda;
% Adimensionalizziamo con grandezze ala
Cm_CG_tot=M_CG_tot/(0.5*rho*(norm(U_Inf)^2)*(c_media1^2)*config.SemiSpan(1)*2);

