function R = raiosModTerra(modTerra,lambda)

    e = modTerra.e;
    R_0 = modTerra.R_0;

    R_E = R_0*(1 + e*(sin(lambda))^2);            % raio leste-oeste
    R_N = R_0*(1 - e*(2 - 3*(sin(lambda))^2));    % raio norte-sul
    R_e = R_0*(1 - e*(sin(lambda))^2);
    
    R = [R_N, R_E, R_e];
end

