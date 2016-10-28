% PM. Practica 2: Implementacio del simplex primal
% Carlos Segarra Gonzalez, Jose A. Ballester Huesca
% Arxiu d'execucio

% Lectura de les dades (A, b, c)
c = dlmread('entrada_deg3/c.dat');
A = dlmread('entrada_deg3/A.dat');
b = dlmread('entrada_deg3/b.dat')';

% Definicio del punter d'escriptura
global fileID;
fileID = fopen('entrada_deg3/sortida2.out', 'w');
fprintf(fileID, '%s\n', 'PM. Practica 2 - Carlos Segarra, Jose Ballester');
fprintf(fileID, '%s\n\n', 'Inici simplex primal amb regla de Bland.');

% Inici fase I
fprintf(fileID, '\t%s\n', 'Fase I');
niterI = 0;
[A, b, VB, VN, z, niterI] = SimplexFaseI(niterI, A, b);
fprintf(fileID, '\t%s', 'Final de la fase I amb ');
fprintf(fileID, '%d', niterI);
fprintf(fileID, '%s\n', ' iteracions.');

% Si l'optim de la fase I es 0 <=> problema factible, fase II
if(abs(z) < 10e-12) % Per questions de precisio
    fprintf(fileID, '\t%s\n', 'Fase II');
    B = A(:, VB);
    invB = inv(B);
    xB = invB*b;
    z = c(VB)*xB;

    % Per defecte, 0: seguim iterant; 1: optim; 2: il.limitat
    res = 0;
    
    niterII = niterI; % Partim del nombre d'iteracions ja fetes
    while(res == 0) 
        niterII = niterII+1;
        fprintf(fileID, '\t\t%s', 'Iteracio: ');
        fprintf(fileID, '%d ', niterII);
        [r, VB, VN, xB, z, res] = ItSimplex (A, VB, VN, xB, z, c);
    end
    
    % Arribem al final
    fprintf(fileID, '\t%s %d %s\n\n', 'Final de la fase II amb', niterII-niterI, 'iteracions.');
    
    % Si hi ha optim
    if(res == 1)
        fprintf(fileID, '%s\n', 'Solucio optima:');
        fprintf(fileID, '\t%s', 'VB = ');
        fprintf(fileID, ' %d ', VB);
        fprintf(fileID, '\n\t%s', 'xB = ');
        fprintf(fileID, '%12.6f', xB);
        fprintf(fileID, '\n\t%s', 'z  = ');
        fprintf(fileID, '%12.6f', z);
        fprintf(fileID, '\n\t%s', 'r  = ');
        fprintf(fileID, '%12.6f', r);
    else
        % Problema il.limitat
        if(res == 2)
        fprintf(fileID, '%s', 'El problema es il.limitat.');
        end
    end
end

fprintf(fileID, '%s\n', '');
fclose(fileID);