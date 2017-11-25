%% Bearbeitungsbogen (-zu bearbeiten-)
clear all
close all
clc

%% Parameter
g = 9.81;                           % m/s^2         %% Erdbeschleunigung
l = 0.2;                            % m             %% Laenge l
AnzahlSchritte = 10000;                             %% Anzahl der Zeitschritte
ts=5;                              % s              %% Laufzeit
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
disp('Berechnungszeit fuer die Loesung der Gleichung mit numerischen Verfahren:');
toc

%% Aufgabe 2: exakte Loesung
x_exakt = [pi/6*cos(sqrt(g/l)*t);-pi/6*sqrt(g/l)*sin(sqrt(g/l)*t)];

tic
%% Aufgabe 6: fehler

l_Eu_expl = x_Eu_expl(2,:)-x_exakt(2,:);
l_Eu_impl = x_Eu_impl(2,:)-x_exakt(2,:);
l_RuKu = x_RuKu(2,:)-x_exakt(2,:);

for i=1:AnzahlSchritte+1
    temp_Eu_ex = x_Eu_expl(1,1:i)-x_exakt(1,1:i);
    temp_Eu_im = x_Eu_impl(1,1:i)-x_exakt(1,1:i);
    temp_RuKu = x_RuKu(1,1:i)-x_exakt(1,1:i);
    if max(temp_Eu_ex) == max(abs(temp_Eu_ex))
        e_Eu_expl(i) = max(temp_Eu_ex);
    else
        e_Eu_expl(i) = -max(abs(temp_Eu_ex));
    end
    if max(temp_Eu_im) == max(abs(temp_Eu_im))
        e_Eu_impl(i) = max(temp_Eu_im);
    else
        e_Eu_impl(i) = -max(abs(temp_Eu_im));
    end
    if max(temp_RuKu) == max(abs(temp_RuKu))
        e_RuKu(i) = max(temp_RuKu);
    else
        e_RuKu(i) = -max(abs(temp_RuKu));
    end
end
disp('Berechnungszeit fuer Fehler: ');
toc

%% Aufgabe 7: Plots
figure(1);
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
%axis([0 max(t) -5 5]);
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
%axis([0 max(t) -5 5]);
grid on;
grid minor;
legend('Euler-exp','Euler-imp','Ru-Ku');
xlabel('t');
ylabel('Fehler');
title('Globaler Fehler - t');

tic
%% h Optimierung
t_count=1;
x_Eu_expl2(:,1) = x_0;
t2(1)=0;
while t2(t_count)<ts
    x_Eu_expl2(:,t_count+1)=x_Eu_expl2(:,t_count)+h.*(SystMatr*x_Eu_expl2(:,t_count));
    h_opt = ((2^p(1)-1)/2^p(1)*F_opt/abs(x_Eu_expl2(1,t_count+1)-pi/6*cos(sqrt(g/l)*(t2(t_count)+h))))^(1/p(1))*h;
    x_Eu_expl2(:,t_count+1)=x_Eu_expl2(:,t_count)+h_opt.*(SystMatr*x_Eu_expl2(:,t_count));
    t2(t_count+1)=t2(t_count)+h_opt;
    t_count=t_count+1;
end
figure(2);
plot(t,x_exakt(1,:),'LineWidth',2);
hold on;
plot(t,x_Eu_expl(1,:),'LineWidth',2);
hold on;
plot(t2,x_Eu_expl2(1,:),'LineWidth',2);
hold on;
axis([0 ts -1 1]);
grid on;
grid minor;
legend('Exakt','Euler-exp','Euler-exp-opt');
xlabel('t');
ylabel('\phi');
title('\phi - t');
disp('Berechnungszeit fuer Loesung der Gleichung mit optimiertem h (nur Euler-expl):');
toc

% for i=1:t_count
%     temp_Eu_ex2 = x_Eu_expl2(1,1:i)-pi/6*cos(sqrt(g/l)*(t2(t_count)));
%     if max(temp_Eu_ex2) == max(abs(temp_Eu_ex2))
%         e_Eu_expl2(i) = max(temp_Eu_ex2);
%     else
%         e_Eu_expl2(i) = -max(abs(temp_Eu_ex2));
%     end
% end
% 
% figure(3);
% plot(t,e_Eu_expl,'LineWidth',2);
% hold on;
% plot(t2,e_Eu_expl2,'LineWidth',2);
% hold on;
% grid on;
% grid minor;
% legend('Euler-exp','Euler-exp-opt');
% xlabel('t');
% ylabel('Fehler');
% title('Globaler Fehler - t');