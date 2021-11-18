%%  Projeto Final - Eletrônica de Potência 2
%
%   Descrição: Simulação do conversor CC-CC Boost Quadrático en malha
%              aberta
%   
%   Degrau de D=0.4 para D=0.51 em t=1 segundos
%% SCRIPT
clc  
close all

%% Dados do Step Time e duração de tempo
tf   = 2;
Dt   = 1e-6;
t    = 0:Dt:tf;

t_degrau = 1;
%% Dados Conversor CC-CC Boost Quadrático
Vin  = 12;
D	 = 0.4;
R    = 5;
C1   = 816.493e-6;
C2   = 1.0202e-3;
L1   = 37e-6;
L2   = 153.029e-6;
fch  = 10e3;

D_final = 0.51;

%% Espaço de Estados para chave aberta e fechada
Aon  = [[     0,       0,      0,     0];
        [     0,       0, (1/L2),     0];
        [     0, (-1/C1),      0,     0];
        [     0,       0,      0,  (-1/(C2*R))]];
    
Bon  = [(1/L1); 0; 0; 0];

Aoff = [[      0,       0, (-1/L1),           0];
        [      0,       0,  (1/L2),     (-1/L2)];
        [ (1/C1), (-1/C1),       0,           0];
        [      0,  (1/C2),       0, (-1/(C2*R))]];

Boff  = [(1/L1); 0; 0; 0];

I     = eye(4);

%% Cálculo das matrizes Mon, Non, Moff e Noff
Mon   = (inv((I - (Dt/2)*Aon)))*(I + (Dt/2)*Aon);
Non   = (inv((I - (Dt/2)*Aon)))*((Dt/2)*Bon);

Moff  = (inv((I - (Dt/2)*Aoff)))*(I + (Dt/2)*Aoff);
Noff  = (inv((I - (Dt/2)*Aoff)))*((Dt/2)*Boff);

%% Simulação Computacional utilizando o método de integração trapezoidal
x = zeros(4,1);
X = zeros(4,length(t));
Vin_t_Dt = 0;
j = 1;

vds = dente_de_serra(fch*t);
for i = t(2:end)
    if D > vds(j)
        x = (Mon*x + Non*(Vin + Vin_t_Dt));
    else
        x = (Moff*x + Noff*(Vin + Vin_t_Dt));
    end
    
    j = j + 1;
    
    if i >= t_degrau && D ~= D_final
         D = D_final;
    end
    X(:,j) = x;       

    Vin_t_Dt = Vin;
end

Il1 = X(1,:);
Il2 = X(2,:);
Vc1 = X(3,:);
Vo  = X(4,:);

figure(1)
plot(t, Il1);
xlabel('{\it t} (sec)');
ylabel('{\it I}_{\rm l1}({\it t})');
title('Corrente no Indutor 1');
grid on;

figure(2)
plot(t, Il2);
xlabel('{\it t} (sec)');
ylabel('{\it I}_{\rm l2}({\it t})');
title('Corrente no Indutor 2');
grid on;

figure(3)
plot(t,Vc1);
xlabel('{\it t} (sec)');
ylabel('{\it V}_{\rm c1}({\it t})');
title('Tensão no Capacitor 1');
grid on;

figure(4)
plot(t, Vo);
xlabel('{\it t} (sec)');
ylabel('{\it V}_{\rm o}({\it t})');
title('Tensão na Saída');
grid on;

%% Extração dos Valores para a aproximação em resposta de segunda ordem
index_t_degrau = find(t == t_degrau);
Vo_max = max(Vo(index_t_degrau:end)); % Tensão de Pico
Tpico     = t(index_t_degrau + find(Vo(index_t_degrau:end) == Vo_max));  % Instante de Pico para o degrau
Vo_antes_degrau = Vo(index_t_degrau); % Tensão em regime permanente ANTES do Degrau
Vo_perm = Vo(end); % Tensão em regime permanente DEPOIS do Degrau

fprintf("\n\nDADOS DA RESPOSTA NA RAZÃO CÍCLICA");
fprintf("\n\nVo(max): %f\n", Vo_max);
fprintf("T(pico): %f\n", Tpico);
fprintf("V(inicial): %f\n", Vo_antes_degrau);
fprintf("T(inicial): %f\n", t_degrau);
fprintf("V(final): %f\n\n", Vo_perm);
