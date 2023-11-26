function [Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo_Input)

NPannelli = length(Corpo_Input.x)-1;

Centro = zeros(NPannelli, 2);
Normale = zeros(NPannelli, 2);
Tangente = zeros(NPannelli, 2);
Estremo_1 = zeros(NPannelli, 2);
Estremo_2 = zeros(NPannelli, 2);
alpha = zeros(NPannelli, 1);
L2G_TransfMatrix = zeros(NPannelli, 2, 2);
G2L_TransfMatrix = zeros(NPannelli, 2, 2);



for i = 1:NPannelli
    
    Centro(i, 1) = ...
    Centro(i, 2) = ...
    
    Estremo_1(i, 1) = ...
    Estremo_1(i, 2) = ...
    
    Estremo_2(i, 1) = ...
    Estremo_2(i, 2) = ...
    
    dy = ...
    dx = ...
    
    % compute angle and matrices
    angle = atan2(dy, dx);
    
    L2G_TransfMatrix(i, :, :) = [cosAngle ,  -sinAngle;
                                 sinAngle,  cosAngle];
                             
    G2L_TransfMatrix(i, :, :) = [cosAngle ,  sinAngle;
                                 -sinAngle,  cosAngle];
                             
    Normale(i, 1) = ...
    Normale(i, 2) = ...
    
    Tangente(i, 1) = ...
    Tangente(i, 2) = ...
    
    lunghezza(i) = norm(Estremo_2(i, :) - Estremo_1(i, :));
    
end


