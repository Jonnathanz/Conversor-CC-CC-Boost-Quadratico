function y = dente_de_serra(t)
%DENTE_DE_SERRA Summary of this function goes here
%   Gera uma fun��o dente de serra com a frequ�ncia de 1 Hz e amplitude
%   unit�ria.
    
    y = 1 - rem(t, 1);
end

