function [q] = normalizacaoIterativa(q0)
% Inicia com DCM computada pelo algoritmo de integração em um dado instante

antes = q0;
n = 1;
while ( (abs(norm(antes) - 1) > 1e-8) && (n < 1e2) )
    depois = (3/2)*antes-(1/2)*quatProd(quatProd(antes,quatConj(antes)),antes);
    antes = depois;
    n = n+1;
end

% Termina com o quaternion unitário q_opt
q = antes;

end

