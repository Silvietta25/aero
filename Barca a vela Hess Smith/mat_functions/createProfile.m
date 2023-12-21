function [x,y] = createProfile(Profilo,NPannelli,Chord,Flap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here




        % Sfruttiamo XFoil per la creazione del profilo, sperando che sia
        % all'interno del database

        fileID = fopen('XFoilInput.txt','w');
        fprintf(fileID, ['naca ' ' '  Profilo, '\n\n']);
        fprintf(fileID,'pane\n\n');

        fprintf(fileID,'gdes\n');

        if nargin==4
           fprintf(fileID,'flap 0.75 0 10\n');
           % nell'ordine sono: coord_x (è sensato sceglierla a 3/4 
           % della corda) e coord_y cerniera flap, angolo di
           % deflessione flap verso il basso in gradi
        end

        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPannelli+1)  '\n\n\n']);

        filename = strcat('NACA_', Profilo, '.dat');

        fprintf(fileID, ['save ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilInput.txt > /dev/null 2>&1");
%         Str2Exec = strcat("xfoil < XFoilInput.txt ");

        system(Str2Exec);

        Corpo = importXfoilProfile(filename);
        
        % Prima flippa i vettori
        x = flipud(Corpo.x);
        y = flipud(Corpo.y);
        
        x = x.*Chord;
        y = y.*Chord;



end