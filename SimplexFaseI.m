% PM. Practica 2: Implementacio del simplex primal
% Carlos Segarra Gonzalez, Jose A. Ballester Huesca
% Funcio de la fase I del simplex primal

function [A3, b3, VB2, VN2, z, niter] = SimplexFaseI (niter, A, b)
    % Adquirim el punter d'escriptura com a variable global
    global fileID;
    
    % Definicio de mides
    m = length(A(:,1));
    n = length(A);
    
    % Comprovem que b sigui mes gran o igual que 0
    % en cas contrari, invertim de signe la fila que tingui b_i negativa
    for i = 1:m
        if(b(i) < 0)
            A(i, :) = -A(i, :);
            b(i) = -b(i);
        end
    end
    
    % Observem quines variables artificials no cal crear
    noCalCrear = [];
    bCobertes = [];
    for i = 1:n 
        bUtil = -1;
        nZeros = 0;
        nUns = 0;
        j = 1;
        while(j <= m && nUns <= 1)
            if(A(j, i) > 0) % Podem usar com a var basica inicial
                nUns = nUns+1;
                bUtil = j;
            end
            if(A(j, i) == 0)
                nZeros = nZeros+1;
            end
            j = j+1;
        end
        if(nZeros == m-1 && nUns == 1) 
            % Aleshores la columna es util
            noCalCrear = [noCalCrear, i];
            bCobertes = [bCobertes, bUtil];
        end
    end
    
    VB = noCalCrear;
    novaVB = n+1;
    A2 = A;
    
    for(i = 1:m)
        if(any(bCobertes == i) == 0) 
            % la fila i-essima de b no esta "coberta" per cap variable natural
            % aleshores, cal crear-ne una d'artificial
            colNova = zeros(length(A2(:, 1)), 1);
            colNova(i, 1) = 1;
            A2 = [A2, colNova];
            VB = [VB, novaVB];
            novaVB = novaVB+1;
        end
    end
    
    % Partim de la SBF trivial per a la fase I
    VN = [];
    c = [];
    % construim VN a partir de les variables que no son basiques
    for i = 1:novaVB-1
        if(any(VB == i) == 0)
            VN = [VN, i];
            c = [c, 0];
        else
            c = [c, 1];
        end
    end

    B = A2(:, VB);
    invB = inv(B);
    xB = invB*b;
    z = c(VB)*xB;
    
    % Per defecte, seguim iterant
    res = 0;
    
    % Iterem fent servir el mateix simplex
    while(res == 0)
        niter = niter+1;
        fprintf(fileID, '\t\t%s', 'Iteracio: ');
        fprintf(fileID, '%d ', niter);
        [r, VB, VN, xB, z, res] = ItSimplex (A2, VB, VN, xB, z, c);
    end
    
    % Mirem la possibilitat que el problema sigui infactible (z* > 0)
    if(abs(z) > 10e-12) % Per questions de precisio
        fprintf(fileID, '\t%s\n', 'El problema es infactible.');
    else
        % Suposem que no es degenerat; ho comprovem
        esDegenerat = 0;
        for i = 1:length(VB)
            if(VB(i) > n)
                % Es degenerat si te alguna variable artificial basica
                esDegenerat = 1;
            end
        end
        
        if(esDegenerat == 1)
            % Si es degenerat, eliminem files redundants o be canviem
            % variables basiques artificials per no artificials
            for l = 1:m
                if(l <= length(VB) && VB(l) > n)
                    fprintf(fileID, '\t%s %d %s\n', 'Problema degenerat: la variable artificial', VB(l), 'es basica despres de la fase I.');
                    B = A2(:, VB);
                    aux = inv(B)*A;
                    
                    % Per defecte, suposem que totes les columnes valen 0 a
                    % la fila l-essima
                    toteszero = 1;
                    for i = 1:length(aux)
                        if(aux(l, i) ~= 0 && any(VB == i) == 0 && toteszero == 1) 
                            % Si ja hem trobat una que valgui zero, podem
                            % substituir la var. artificial VB(l) per i
                            toteszero = 0;
                            fprintf(fileID, '\t%s %d %s %d%s\n', 'Hem substituit la variable artificial', VB(l), 'per la no artificial', i, '.');
                            VB(l) = i;
                        end
                    end
                    if(toteszero == 1)
                        % Si totes son zero, vol dir que la fila es
                        % redundant i per tant pot ser eliminada
                        A = [A(1:l-1, :); A(l+1:m, :)];
                        b = [b(1:l-1); b(l+1:m)];
                        VB = [VB(1:l-1), VB(l+1:m)];
                        fprintf(fileID, '\t%s %d %s\n', 'Hem eliminat la constriccio num.', l, 'per ser redundant.');
                    end
                end
            end
            % Fi del cas degenerat
        end
    end
    
    % Actualitzacions de cara a la fase II
    VB2 = VB;
    VN2 = [];
    for i = 1:n
        if(any(VB == i) == 0)
            VN2 = [VN2, i];
        end
    end
    A3 = A;
    b3 = b;
end