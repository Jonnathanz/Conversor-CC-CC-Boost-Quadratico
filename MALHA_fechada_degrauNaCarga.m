%%  Projeto Final - Eletrônica de Potência 2
%
%   Descrição: Simulação do conversor CC-CC Boost Quadrático en malha
%              fechada
%   
%   Degrau de D=0.4 para D=0.51 em t=0.2 segundos
%% SCRIPT
clc  
close all
clear all

%% Dados do Step Time e duração de tempo
tf   = 2;
Dt   = 1e-6;
t    = 0:Dt:tf;

t_degrau = 1;
R_novo = 5;
%% Dados Conversor CC-CC Boost Quadrático
Vin  = 12;
D	 = 0.5;
R    = 10;
C1   = 816.493e-6;
C2   = 1.0202e-3;
L1   = 37e-6;
L2   = 153.029e-6;
fch  = 10e3;

%% Dados do controlador PID
kd = 1*4.66789e-7;
kp = 3.46173e-5;
ki = 0.28858;

Vref = 50;
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
d = zeros(1,length(t));
Vin_t_Dt = 0;

j = 1;

e_antigo = 0;
e_atual = 0;

P = 0;
Itg = 0;
De= 0;

for i = t(2:end)
    j = j + 1;
    
    e_atual = Vref - x(4);
    
    P  = kp*e_atual;
    Itg  = Itg + ki*Dt*(e_atual + e_antigo)/2;
    De = kd*(e_atual - e_antigo)/Dt;
    
    D = P + Itg + De;
    
    if D >= 0.7
        D = 0.7;
    elseif D < 0
         D = 0;
    end
    
    d(j) = D;
    
    if D >= dente_de_serra(fch*i)
        x = (Mon*x + Non*(Vin + Vin_t_Dt));
    else
        x = (Moff*x + Noff*(Vin + Vin_t_Dt));
    end
    
    if i > t_degrau && R ~= R_novo
        R     = R_novo;
        
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
        
        Mon   = (inv((I - (Dt/2)*Aon)))*(I + (Dt/2)*Aon);
        Non   = (inv((I - (Dt/2)*Aon)))*((Dt/2)*Bon);
        Moff  = (inv((I - (Dt/2)*Aoff)))*(I + (Dt/2)*Aoff);
        Noff  = (inv((I - (Dt/2)*Aoff)))*((Dt/2)*Boff);
    end

    X(:,j) = x;     
    
    Vin_t_Dt = Vin;
    e_antigo = e_atual;
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
xlim([0.3, tf]);
% ylim([0 60]);
grid on;

figure(2)
plot(t, Il2);
xlabel('{\it t} (sec)');
ylabel('{\it I}_{\rm l2}({\it t})');
title('Corrente no Indutor 2');
xlim([0.3, tf]);
ylim([0 60]);
grid on;

figure(3)
plot(t,Vc1);
xlabel('{\it t} (sec)');
ylabel('{\it V}_{\rm c1}({\it t})');
title('Tensão no Capacitor 1');
xlim([0.3, tf]);
ylim([0 60]);
grid on;

figure(4)
plot(t, Vo);
xlabel('{\it t} (sec)');
ylabel('{\it V}_{\rm o}({\it t})');
title('Tensão na Saída');
xlim([0.3, tf]);
ylim([0 60]);
grid on;

figure(5)
plot(t, d);
xlabel('{\it t} (sec)');
ylabel('D');
title('Razão Cíclica');
xlim([0.3, tf]);
ylim([0 0.7]);
grid on;

% min(Vo(900000:end))