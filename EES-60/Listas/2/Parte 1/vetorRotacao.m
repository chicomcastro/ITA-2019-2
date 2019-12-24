function [desalinhamento] = vetorRotacao(q_A,q_B)
    q = quatProd(q_B, quatInv(q_A));
    
    Psi = quatToRot(q);
    desalinhamento = Psi;
end

