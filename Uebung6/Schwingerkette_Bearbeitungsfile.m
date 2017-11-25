%% Bearbeitungsbogen (-freiwillig-)
clear all;
close all;
clc;

%% Parameter
m=0.5;                      % kg
c=500;                      % N/m
AnzahlSchritte=10000;
ts=1;                       %s
h=ts/AnzahlSchritte;        %s

%% Systemmatrix
% m*x_1" + 2k*x_1 - k*x_2 = 0
% m*x_2" - k*x_1 + 2k*x_1 = 0
A_Sys=[0 0 1 0; 0 0 0 1; -2*c/m c/m 0 0; c/m -2*c/m 0 0];

%% Anfangsbedingungen            
x_1(1)=0.01;
x_2(1)=0.01;
xa_1(1)=0;
xa_2(1)=0;
X(:,1)=[x_1(1);x_2(1);xa_1(1);xa_2(1)];
t(1)=0;

%% Aufgabe 2: Eigenwerte
[Eigenvektoren, Eigenwerte] = eig(A_Sys);

%% Aufgabe 3: Verfahren von Heun            
for n = 1:1:AnzahlSchritte
    k1(:,n) = h.*(A_Sys*X(:,n));
    k2(:,n) = h.*(A_Sys*(X(:,n)+k1(:,n)));
    X(:,n+1) = X(:,n)+1/2.*(k1(:,n)+k2(:,n));
    
    t(n+1)=t(n)+h;
end
 
%% Aufgabe 4: Berechnung der Eigenfrequenz bei zunehmdener Federsteifigkeit
figure(1);
subplot('position',[0.05,0.55,0.4,0.4]);
plot(t,X(1,:),'LineWidth',2);
hold on;
plot(t,X(2,:),'LineWidth',2);
hold on;
grid on;
grid minor;
xlabel('t');
ylabel('x');
title('\itt - x');
legend('x1','x2');
axis([0 ts -0.012 0.012]);

subplot('position',[0.55,0.55,0.4,0.4]);
plot(t,X(3,:),'LineWidth',2);
hold on;
plot(t,X(4,:),'LineWidth',2);
hold on;
grid on;
grid minor;
xlabel('t');
ylabel('v');
title('\itt - v');
legend('v1','v2');
%axis([0 ts -0.4 0.4]);

c2=1:1:1000;
for i=1:1:1000
    A_Sys2=[0 0 1 0; 0 0 0 1; -2*i/m i/m 0 0; i/m -2*i/m 0 0];
    lamda=eig(A_Sys2);
    EigenFre1(i)=abs(lamda(1));
    EigenFre2(i)=abs(lamda(3));
end

subplot('position',[0.05,0.05,0.9,0.4]);
plot(c2,EigenFre1,'LineWidth',2);
hold on;
plot(c2,EigenFre2,'LineWidth',2);
hold on;
grid on;
grid minor;
xlabel('Federsteifigfkeit');
ylabel('Eigenfrequenze');
title('\itSteifigkeit - Eigenfrequenze');
legend('Eigenfrequenze 1','Eigenfrequenze 2');