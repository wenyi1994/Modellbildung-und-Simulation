%% MathPendel_Bearbeitungsfile.m
%  - Systemmatrix erstellen
%  - Anfangsbedingungen einsetzen
%  - Loesung der Gleichung mit verschiedenen numerischen Verfahren
%    berechnen
%  - Lakale und Globale Fehler berechnen
%  - Zeitschritte optimieren und vergleichen
%  ----------------------------------------------------------------------
%  Output:
%       - Berechnungszeit fuer die Loesung der Gleichung(en) mit
%         numerischen Verfahren
%       - Berechnungszeit fuer Fehler
%       - Berechnungszeti fuer die Loesung der Gleichung(en) mit
%         optimierten h (Zeitschrittweite)
%       - Globale Diskretisierungsfehler nach max. Anzahl der Schritte
%       - Figure 1: Numerische Verfahren mit fixen 'h'
%           + Zeitlicher Verlauf von Winkel
%           + Zeitlicher Verlauf von Winkelgeschwindigkeit
%           + Zeitlicher Verlauf von lokalen Fehler
%           + Zeitlicher Verlauf von glabalen Fehler
%       - Figure 2: Numerische Verfahren mit optimierten 'h'
%           + Zeitlicher Verlauf von Winkel mit Euler-Explizit-Verfahren
%           + Zeitlicher Verlauf von Winkel mit Euler-Implizit-Verfahren
%           + Zeitlicher Verlauf von Winkel mit Runge-Kutta-Verfahren
%           + Zeitlicher Verlauf von Winkel
%  ----------------------------------------------------------------------
%  Wen Yi, Karlsruhe Institut of Technology
%  yi.wen@student.kit.edu
%  github.com/wenyi1994/Modellbildung-und-Simulation
%  2017/11/28

%% Bearbeitungsbogen (-zu bearbeiten-)
clear all
close all
clc

%% Parameter
g = 9.81;                           % m/s^2         %% Erdbeschleunigung
l = 0.2;                            % m             %% Laenge l
AnzahlSchritte = 10000;                             %% Anzahl der Zeitschritte
ts=5;                               % s             %% Laufzeit
F_opt=0.0001;                                       %% vorgegebener Fehler
h = ts/AnzahlSchritte;                              %% Zeitschrittweite
x_0 = [pi/6;0];                                     %% Anfangsbedingungen in Zustandsvektor
t(1) = 0;                           % s             %% Initialisierung der Zeit
p=[1;4];                                            %% Fehlerordnung

x_Eu_expl(:,1) = x_0;                               %% Anfangsbedingungen fuer Euler explizit
x_Eu_impl(:,1) = x_0;                               %% Anfangsbedingungen fuer Euler implizit
x_RuKu(:,1) = x_0;                                  %% Anfangsbedingungen fuer Runge-Kutta-Verfahren

SystMatr=[0, 1; -g./l, 0];                          %% Systemmatrix
  

tic
for n = 1:AnzahlSchritte
    %% Aufgabe 3: Euler explizit
    x_Eu_expl(:,n+1)=x_Eu_expl(:,n)+h.*(SystMatr*x_Eu_expl(:,n));
       
    %% Aufgabe 4: Euler implizit
    x_Eu_impl(:,n+1) =(eye(2)-h.*SystMatr)\x_Eu_impl(:,n);
    
    %% Aufgabe 5: Runge-Kutta-Verfahren    
    k1(:,n)=h.*(SystMatr*x_RuKu(:,n));
    k2(:,n)=h.*(SystMatr*(x_RuKu(:,n)+0.5.*k1(:,n)));
    k3(:,n)=h.*(SystMatr*(x_RuKu(:,n)+0.5.*k2(:,n)));
    k4(:,n)=h.*(SystMatr*(x_RuKu(:,n)+k3(:,n)));
    x_RuKu(:,n+1) = x_RuKu(:,n)+1/6.*(k1(:,n)+2*k2(:,n)+2*k3(:,n)+k4(:,n));
    
    t(n+1)=t(n)+h;
end
disp('Berechnungszeit: Numerische Verfahren mit fixen ''h''');
disp(['                                    ',num2str(toc),' seconds']);
disp('----------------------------------------------------------------------');

%% Aufgabe 2: exakte Loesung
x_exakt = [pi/6*cos(sqrt(g/l)*t);-pi/6*sqrt(g/l)*sin(sqrt(g/l)*t)];

tic
%% Aufgabe 6: fehler
l_Eu_expl(1)=0;
l_Eu_impl(1)=0;
l_RuKu(1)=0;
for i=1:AnzahlSchritte
    l_Eu_expl(i+1) = (x_Eu_expl(1,i+1)-x_Eu_expl(1,i))/h - (x_exakt(1,i+1)-x_exakt(1,i))/h;
    l_Eu_impl(i+1) = (x_Eu_impl(1,i+1)-x_Eu_impl(1,i))/h - (x_exakt(1,i+1)-x_exakt(1,i))/h;
    l_RuKu(i+1) = (x_RuKu(1,i+1)-x_RuKu(1,i))/h - (x_exakt(1,i+1)-x_exakt(1,i))/h;
end

for i=1:AnzahlSchritte+1
    temp_Eu_ex = x_Eu_expl(1,1:i)-x_exakt(1,1:i);
    temp_Eu_im = x_Eu_impl(1,1:i)-x_exakt(1,1:i);
    temp_RuKu = x_RuKu(1,1:i)-x_exakt(1,1:i);
    
    temp_1 = max(temp_Eu_ex);
    temp_2 = max(abs(temp_Eu_ex));
    e_Eu_expl(i) = [temp_1 == temp_2, temp_1 ~= temp_2] * [temp_2; -temp_2];
    temp_1 = max(temp_Eu_im);
    temp_2 = max(abs(temp_Eu_im));
    e_Eu_impl(i) = [temp_1 == temp_2, temp_1 ~= temp_2] * [temp_2; -temp_2];
    temp_1 = max(temp_RuKu);
    temp_2 = max(abs(temp_RuKu));
    e_RuKu(i) = [temp_1 == temp_2, temp_1 ~= temp_2] * [temp_2; -temp_2];
end
disp('Berechnungszeit: Fehler');
disp(['                                    ',num2str(toc),' seconds']);
disp('----------------------------------------------------------------------');
disp(['Globale Diskretisierungsfehler nach ',num2str(AnzahlSchritte),' Schritte:']);
disp(' ');
disp(['                Euler-Explizit-Verfahren        |       ',num2str(e_Eu_expl(AnzahlSchritte+1))]);
disp(['                Euler-Implizit-Verfahren        |       ',num2str(e_Eu_impl(AnzahlSchritte+1))]);
disp(['                  Runge-Kutta-Verfahren         |       ',num2str(e_RuKu(AnzahlSchritte+1))]);
disp('----------------------------------------------------------------------');

%% Aufgabe 7: Plots
figure('Name','Numerische Verfahren mit fixen ''h''','NumberTitle','off');
subplot('position',[0.05,0.55,0.4,0.4]);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t,x_Eu_expl(1,:),'LineWidth',2);
hold on;
plot(t,x_Eu_impl(1,:),'LineWidth',2);
hold on;
plot(t,x_RuKu(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Euler-exp','Euler-imp','Ru-Ku');
xlabel('t');
ylabel('\phi');
title('\phi - t');

subplot('position',[0.55,0.55,0.4,0.4]);
plot(t,x_exakt(2,:),'LineWidth',2);
hold on;
plot(t,x_Eu_expl(2,:),'LineWidth',2);
hold on;
plot(t,x_Eu_impl(2,:),'LineWidth',2);
hold on;
plot(t,x_RuKu(2,:),'LineWidth',2);
hold on;
axis([0 ts -5 5]);
grid on;
grid minor;
legend('Exakt','Euler-exp','Euler-imp','Ru-Ku');
xlabel('t');
ylabel('\phi''');
title('\phi'' - t');

%% Aufgabe 8: Plots
subplot('position',[0.05,0.05,0.4,0.4]);
plot(t,l_Eu_expl,'LineWidth',2);
hold on;
plot(t,l_Eu_impl,'LineWidth',2);
hold on;
plot(t,l_RuKu,'LineWidth',2);
hold on;
grid on;
grid minor;
legend('Euler-exp','Euler-imp','Ru-Ku');
xlabel('t');
ylabel('Fehler');
title('Lokaler Fehler - t');

subplot('position',[0.55,0.05,0.4,0.4]);
plot(t,e_Eu_expl,'LineWidth',2);
hold on;
plot(t,e_Eu_impl,'LineWidth',2);
hold on;
plot(t,e_RuKu,'LineWidth',2);
hold on;
grid on;
grid minor;
legend('Euler-exp','Euler-imp','Ru-Ku');
xlabel('t');
ylabel('Fehler');
title('Globaler Fehler - t');

tic
%% h Optimierung
x_Eu_expl_h(:,1) = x_0;
x_Eu_impl_h(:,1) = x_0;
x_RuKu_h(:,1) = x_0;

t_count = 1;
t2(1) = 0;
t3(1) = 0;
t4(1) = 0;
while t2(t_count)<ts
    %% Euler-Exp mit optimierten h
    x_Eu_expl_h(:,t_count+1)=x_Eu_expl_h(:,t_count)+h.*(SystMatr*x_Eu_expl_h(:,t_count));
    h_opt=h;
    delta=abs(x_Eu_expl_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t2(t_count)+h_opt)));
    if delta > F_opt
        h_opt = ((2^p(1)-1)/2^p(1)*F_opt/delta)^(1/p(1))*h_opt;
        x_Eu_expl_h(:,t_count+1)=x_Eu_expl_h(:,t_count)+h_opt.*(SystMatr*x_Eu_expl_h(:,t_count));

        delta=abs(x_Eu_expl_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t2(t_count)+h_opt)));
    end
    t2(t_count+1)=t2(t_count)+h_opt;
    t_count=t_count+1;
end
t_count = 1;
while t3(t_count)<ts
    %% Euler-Imp mit optimierten h
    x_Eu_impl_h(:,t_count+1) =(eye(2)-h.*SystMatr)\x_Eu_impl_h(:,t_count);
    h_opt=h;
    delta=abs(x_Eu_impl_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t3(t_count)+h_opt)));
    if delta > F_opt
        h_opt = ((2^p(1)-1)/2^p(1)*F_opt/delta)^(1/p(1))*h_opt;
        x_Eu_impl_h(:,t_count+1) =(eye(2)-h_opt.*SystMatr)\x_Eu_impl_h(:,t_count);

        delta=abs(x_Eu_impl_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t3(t_count)+h_opt)));
    end
    t3(t_count+1)=t3(t_count)+h_opt;
    t_count=t_count+1;
end
t_count = 1;
while t4(t_count)<ts
    %% Runge-Kutta mit optimierten h
    k1h(:,t_count)=h.*(SystMatr*x_RuKu_h(:,t_count));
    k2h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+0.5.*k1h(:,t_count)));
    k3h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+0.5.*k2h(:,t_count)));
    k4h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+k3h(:,t_count)));
    x_RuKu_h(:,t_count+1) = x_RuKu_h(:,t_count)+1/6.*(k1h(:,t_count)+2*k2h(:,t_count)+2*k3h(:,t_count)+k4h(:,t_count));
    h_opt=h;
    delta=abs(x_RuKu_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t4(t_count)+h_opt)));
    if delta > F_opt
        h_opt = ((2^p(1)-1)/2^p(1)*F_opt/delta)^(1/p(1))*h_opt;
        k1h(:,t_count)=h.*(SystMatr*x_RuKu_h(:,t_count));
        k2h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+0.5.*k1h(:,t_count)));
        k3h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+0.5.*k2h(:,t_count)));
        k4h(:,t_count)=h.*(SystMatr*(x_RuKu_h(:,t_count)+k3h(:,t_count)));
        x_RuKu_h(:,t_count+1) = x_RuKu_h(:,t_count)+1/6.*(k1h(:,t_count)+2*k2h(:,t_count)+2*k3h(:,t_count)+k4h(:,t_count));

        delta=abs(x_RuKu_h(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t4(t_count)+h_opt)));
    end
    t4(t_count+1)=t4(t_count)+h_opt;
    t_count=t_count+1;
end
disp('Berechnungszeit: Numerische Verfahren mit optimierten ''h''');
disp(['                                    ',num2str(toc),' seconds']);

%% Plots mit optimierten h
figure('Name','Numerische Verfahren mit optimierten ''h''','NumberTitle','off');
subplot('position',[0.05,0.55,0.4,0.4]);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t,x_Eu_expl(1,:),'LineWidth',2);
hold on;
plot(t2,x_Eu_expl_h(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Euler-exp','Euler-exp-opt');
xlabel('t');
ylabel('\phi');
title({'Euler-Explizit';'\phi - t'});

subplot('position',[0.55,0.55,0.4,0.4]);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t,x_Eu_impl(1,:),'LineWidth',2);
hold on;
plot(t3,x_Eu_impl_h(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Euler-imp','Euler-imp-opt');
xlabel('t');
ylabel('\phi');
title({'Euler-Implizit';'\phi - t'});

subplot('position',[0.05,0.05,0.4,0.4]);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t,x_RuKu(1,:),'LineWidth',2);
hold on;
plot(t4,x_RuKu_h(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Runge-Kutta','Runge-Kutta-opt');
xlabel('t');
ylabel('\phi');
title({'Runge-Kutta';'\phi - t'});

subplot('position',[0.55,0.05,0.4,0.4]);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t2,x_Eu_expl_h(1,:),'LineWidth',2);
hold on;
plot(t3,x_Eu_impl_h(1,:),'LineWidth',2);
hold on;
plot(t4,x_RuKu_h(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Euler-exp','Euler-imp','Runge-Kutta');
xlabel('t');
ylabel('\phi');
title({'Numerische Verfahren mit optimierten ''h''';'\phi - t'});
