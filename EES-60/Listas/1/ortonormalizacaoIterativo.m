function [D] = ortonormalizacaoIterativo(D0)
% Inicia com DCM computada pelo algoritmo de integração em um dado instante

antes = D0;
n = 1;
while ((abs(det(antes) - 1) > 1e-8) && (n < 1e2))
    depois = (3/2)*antes-(1/2)*antes*antes'*antes;
    antes = depois;
    n = n+1;
end

% Termina com DCM sendo a matriz ortonormal Dopt que minimiza a norma euclidiana
% entre D0 e Dopt
D = antes;

end

