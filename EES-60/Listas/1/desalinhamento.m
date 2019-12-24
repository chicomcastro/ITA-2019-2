function [desalinhamento, vetorDesalinhamento] = desalinhamento(A,B)
    D = A'*B;
    E = D - eye(length(D));
    E_par = (E+E')/2;
    E_impar = (E-E')/2;
    
    Psi = [
        ( abs(E_impar(2,3)) + abs(E_impar(3,2)) ) / 2;
        ( abs(E_impar(1,3)) + abs(E_impar(3,1)) ) / 2;
        ( abs(E_impar(1,2)) + abs(E_impar(2,1)) ) / 2
        ];
    vetorDesalinhamento = Psi;
    desalinhamento = norm(Psi);
end

