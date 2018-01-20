function [u,fa] = fem_Bearbeitungsfile(M)
% fem_Bearbeitungsfile.m
% Modellbildung und Simulation
% Übung zum Kapitel: Numerische Lösungsverfahren für partielle DGLen: Finite-Elemente-Methode
% ----------------------------------------------------------------------
% [u,fa] = fem_Bearbeitungsfile(M)
% ----------------------------------------------------------------------
% input:
% M                     Anzahl der Elemente
% ----------------------------------------------------------------------
% output:
% u                     Verschiebung verschiedener Stelle (x)
% fa                    Reaktionskraft an Stabrand
% ----------------------------------------------------------------------
% Wen Yi, Karlsruhe Institut of Technology
% yi.wen@student.kit.edu
% 2018/01/20
close all;

%% Parameter
L = 1;
F = 20000;
E = 70e9;
AnzahlElemente = M;
LElemente = L/AnzahlElemente;   %% Länge eines Elemntes
Knoten = 0:LElemente:L;  %% Knotenvektor
AElement = (10-9*Knoten) * 1e-4;    %% Querschnittsverlauf

%% Aufgabe 1: DGL
%  du/dx = 1/(3500 - 3150*x)---------m
%  du/dx = 1/(3500 - 31.5*x)---------cm

%% Aufgabe 2: Analytische Lösung
    % syms u;
    u = dsolve('Du=1/(3500 - 31.5 * x)','u(0)=0','x')           % Analytische Lösung
    % x = 0:1:100;
    % u = eval(u);
    % figure('name','Finite Element Method','numbertitle','off');
    % plot(x,u,'linewidth',2);
    % hold on;

[x,u] = ode45(@(x,u)1/(3500 - 3150 * x),[0,1],0);           % Numerische Lösung
plot(x,u,'linewidth',2);
hold on;

% %% Aufgabe 3: Formfunktion N
% N = [1 - Knoten'/L, Knoten'/L];
x_zahl = 100;
x = 0:LElemente/x_zahl:LElemente;
N = [1 - x'/LElemente, x'/LElemente];

%% Aufgabe 4: B und D
% B
A = 0.5 * AElement(1:AnzahlElemente) + 0.5 * AElement(2:AnzahlElemente+1);
B = [-1/LElemente, 1/LElemente];

% D
D = E;

%% Aufgabe 5:
for b = 1:1:AnzahlElemente           %% Schleife für die Anzahl an gewählten Elementen
   K{b} = A(b) * LElemente * B' * D * B;
end

%% Aufgabe 6
%% KGesamt
KGesamt = zeros(length(Knoten),length(Knoten));
for b = 1:1:AnzahlElemente
    KGesamt(b:b+1,b:b+1) = KGesamt(b:b+1,b:b+1)+K{b};
end

%% Aufgabe 7
ua = 0;
fb = [zeros(AnzahlElemente-1,1);F];

Kaa = KGesamt(1,1);
Kab = KGesamt(1,2:AnzahlElemente+1);
Kba = KGesamt(2:AnzahlElemente+1,1);
Kbb = KGesamt(2:AnzahlElemente+1,2:AnzahlElemente+1);

ub = Kbb \ fb;
fa = Kaa * ua + Kab * ub;

% %% Aufgabe 8
u = [ua;ub];
stairs(Knoten,u,'linewidth',2);

grid on;
grid minor;
xlabel('x/m');
ylabel('u/m');
legend('Matlab-Numerische-Loesung','FEM-Loesung');


end