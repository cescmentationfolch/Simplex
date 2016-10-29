% simplex fase 1 prova crec que no va

function [VB2, z] = SimplexFaseI (A, b)
    m = length(A(:,1));
    n = length(A);
    Id = zeros(m);
    A = [A, Id];
    VB = n+1:n+m;
    VN = 1:n;
    c = [zeros(1, n), ones(1, m)];
    B = A(:, VB);
    invB = inv(B);
    xB = invB*b;
    z = c(VB)*xB;
    res = 0;
    while(res == 0)
        [VB, VN, xB, z, res] = ItSimplexFaseII (A, VB, VN, xB, z, c);
    end
    z
end
