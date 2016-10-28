% PM. Practica 2: Implementacio del simplex primal
% Carlos Segarra Gonzalez, Jose A. Ballester Huesca
% Iteracio del simplex primal

function [r, VB2, VN2, xB2, z2, res] = ItSimplex (A, VB, VN, xB, z, c)
    global fileID;
    % res es l'indicador del resultat de la iteracio:
    % 0 si ha de seguir iterant, 1 si la solucio es optima, 2 si el
    % problema es il.limitat
    res = 0;
    
    % Definim costos, matrius, costos reduits
    cB = c(1, VB);
    cN = c(1, VN);
    B = A(:, VB);
    AN = A(:, VN);
    invB = inv(B);
    r = cN - cB*invB*AN;
    
    % q es la variable sortint (index absolut)
    % w es la variable sortint (index relatiu dins N)
    % w ens sera util per a substituir la VNB
    q = -1;
    w = -1;
    
    for i = 1:length(r)
        % Si trobem alguna component negativa, seleccionem l'index corresponent
        % si aquest es mes petit que l'agafat actualment (regla de Bland)
        if(r(i) < 0)
            qact = VN(i);
            if(q == -1 || qact < q) 
                w = i; % w = N^{-1}(q)
                q = qact;
            end
        end
    end
    if(q == -1)
        % Aleshores no hem trobat cap component negativa, r >= 0, optim
        res = 1;
        fprintf(fileID, '\t%s', 'res: ');
        fprintf(fileID, '%d', res);
        VB2 = VB;
        VN2 = VN;
        xB2 = xB;
        z2 = z;
    else
        % Si no es optim
        dB = -invB*A(:, q);
        theta = -1;
        p = -1;
        % Anem definint theta i prenem el minim en cada moment
        for i = 1:length(B)
            if(dB(i) < 0)
                thetaact = -xB(i)/dB(i);
                if(theta == -1 || thetaact <= theta) 
                    if(p == -1 || thetaact < theta || VB(i) < VB(p))
                        p = i;
                        % Regla de Bland per a la variable sortint
                    end
                    theta = thetaact;
                end
            end
        end
        if(theta == -1)
            % Problema il.limitat
            res = 2;
            fprintf(fileID, '\t%s', 'res: ');
            fprintf(fileID, '%d', res);
            VB2 = VB;
            VN2 = VN;
            xB2 = xB;
            z2 = z;
        else
            % El problema no es il.limitat
            fprintf(fileID, '\t%s', 'res: ');
            fprintf(fileID, '%d ', res);
            fprintf(fileID, '\t%s', 'q: ');
            fprintf(fileID, '%d ', q);
            fprintf(fileID, '\t%s', 'B(p): ');
            fprintf(fileID, '%d ', VB(p));
            fprintf(fileID, '\t%s', 'theta: ');
            fprintf(fileID, '%12.6f ', theta);
        
            % Actualitzacions
            xB = xB + theta*dB;
            novaVN = VB(p);
            VB(p) = q;
            VN(w) = novaVN;
            xB(p) = theta;
            z = z + theta*r(w);
            fprintf(fileID, '\t%s', 'z: ');
            fprintf(fileID, '%12.6f', z);
            
            z2 = z;
            xB2 = xB;
            VB2 = VB;
            VN2 = VN;
        end
    end
    fprintf(fileID, '%s\n', '');
end