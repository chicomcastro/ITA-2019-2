
deramCerto = 0;
for i = 1:1e2
q = [rand rand rand rand]';
q = q/norm(q)*1.00001;

q_norm = normalizacaoIterativa(q);
if (abs(norm(q_norm) - 1) < 1e-4)
    deramCerto = deramCerto + 1;
end

end

deramCerto