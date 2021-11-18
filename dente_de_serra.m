function y = dente_de_serra(t)
%DENTE_DE_SERRA Summary of this function goes here
%   Gera uma função dente de serra com a frequência de 1 Hz e amplitude
%   unitária.
    
    y = 1 - rem(t, 1);
end

