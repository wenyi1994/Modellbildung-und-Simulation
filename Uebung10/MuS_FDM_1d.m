function [phi_uds_t,phi_cds_t] = MuS_FDM_1d(varargin)
% MuS_FDM_1d.m
% Modellbildung und Simulation
% Übung zum Kapitel: Zeitkontinuierliche Modelle mit verteilten Parametern
% Dr.-Ing. Balazs Pritz, pritz@kit.edu
% 1D Testfall: Konvektions-Diffusionsgleichung mit Dirichlet-Randbedingungen
% Diskretisierungsschema: Finite Differenzen Methode
% Quelle: Joel H. Ferziger, Milovan Peric:
% Numerische Strömungsmechanik, Springer, 2008
% ----------------------------------------------------------------------
% [phi_uds,phi_cds] = MuS_FDM_1d('dt',dt,'nt',nt,'u0',u0,'nx',nx,'rho',rho,'gamma',gamma,'output',n_output,'tol',tol,'init','line')
% [phi_uds,phi_cds] = MuS_FDM_1d('default')
% [phi_uds,phi_cds] = MuS_FDM_1d('default','dt',dt,'nt',nt,...)
% ----------------------------------------------------------------------
% input:
% dt                    Zeitschritt
% maxnt                 Maximale Anzahl von Zeitschritten
%                       (Simulationszeit=timestep*maxnt)
%                       (Es ist auch möglich, dass die Simulation nicht nur
%                       durch maxnt beendet wird, sondern die Änderung der
%                       Lösung von zwei Zeitschritten berechnet wird und
%                       anhand dieser ein Abbruchkriterium definiert wird.)
% u0                    Konvektionsgeschwindigkeit
% nx                    Anzahl der Knotenpunte
% rho                   Dichte, Anzahl oder 'water', 'air', ...
% gamma                 Diffusionskoeffizient für phi
% n_output              Anzahl der insgesamt gegebenen Grafiken
% tol                   Toleranz der Convergenz
% init                  Modul der Initialisierung von phi, 'zero', 'line', ...
% ----------------------------------------------------------------------
% output:
% phi_uds               Lösung der PDGL mit UDS-Verfahren
% phi_cds               Lösung der PDGL mit CDS-Verfahren
% ----------------------------------------------------------------------
% Wen Yi, Karlsruhe Institut of Technology
% yi.wen@student.kit.edu
% 2017/12/24

% Parametern
if strcmp(varargin{1},'default')
    % Default Setting
    dt = 0.01;
    maxnt = 200;
    u0 = 1;
    nx = 11;
    rho = 1.25;
    gamma = 0.01;
    n_output = 50;
    tol = 1e-5;
    init_mod = 'zero';
end
for i = 1:nargin-1
    switch varargin{i}
        case 'dt'
        dt = varargin{i+1};
        case 'nt'
        maxnt = varargin{i+1};
        case 'u0'
        u0 = varargin{i+1};
        case 'nx'
        nx = varargin{i+1};
        case 'rho'
        switch varargin{i+1}
            case 'air'
            rho = 1.25;
            case 'water'
            rho = 998;
            otherwise
            rho = varargin{i+1};
        end
        case 'gamma'
        gamma = varargin{i+1};
        case 'output'
        n_output = varargin{i+1};
        case 'tol'
        tol = varargin{i+1};
        case 'init'
        init_mod = varargin{i+1};
    end
end

nxm1=nx-1;

% Zellengröße
dx=1/(nxm1);

% Initialisierung
phi_uds=zeros(nx,1);
phi_cds=zeros(nx,1);

% Convergency
temp_uds=zeros(nx,1);
temp_cds=zeros(nx,1);
conv_uds = 0;
conv_zeit_uds = 0;
conv_cds = 0;
conv_zeit_cds = 0;

% Randbedingungen für phi
phi_uds(1)=0;
phi_uds(nx)=1;

phi_cds(1)=0;
phi_cds(nx)=1;

switch init_mod
case 'zero'
case 'line'
    for i = 2:nxm1
        phi_uds(i) = (nx - i) / nx * phi_uds(1) + i / nx * phi_uds(nx);
        phi_cds(i) = (nx - i) / nx * phi_cds(1) + i / nx * phi_cds(nx);
    end
case 'random'
    for i = 2:nxm1
        phi_uds(i) = phi_uds(1) * rand + (phi_uds(nx) - phi_uds(1)) * rand;
        phi_cds(i) = phi_cds(1) * rand + (phi_cds(nx) - phi_cds(1)) * rand;
    end
case 'boundaryvalue1'
    for i = 2:nxm1
        phi_uds(i) = phi_uds(1);
        phi_cds(i) = phi_cds(1);
    end
case 'boundaryvalue2'
    for i = 2:nxm1
        phi_uds(i) = phi_uds(nx);
        phi_cds(i) = phi_cds(nx);
    end
otherwise
    warning('Unexpected initialization, ''zero'' is used');
end

% Peclet-Zahl
% Für die Änderung der Pe-Zahl können Sie natürlich beliebig rho, u oder gamma ändern
peclet=rho*u0/gamma;

n=floor(maxnt/n_output);
phi_uds_t = zeros(n_output+1,nx);
phi_cds_t = zeros(n_output+1,nx);


% Die exakte Lösung (wird mit besserer Auflösung gerechnet)
x_n = floor(1/dt);
xel=zeros(x_n+1,1);
el=zeros(x_n+1,1);
for i=1:x_n
    xel(i)=dt*i-dt;
    el(i)=(exp(xel(i)*peclet) - 1.)/(exp(peclet) - 1.);
end
xel(x_n+1) = 1;
el(x_n+1) = 1;

% Zur Beurteilung der Stabilität wird DCFL gerechnet.
% (Wenn Sie eine lokale Verfeinerung verwenden, müssen Sie nach dem größten Wert suchen.)
dcfl= 2*gamma*dt/(dx*dx*rho) + u0*dt/dx;
cfl = u0*dt/dx;
d = gamma*dt/(dx*dx*rho);
t = rho * dx * dx / 2 / gamma;

% Die Position der Stützstellen wird berechnet
x=zeros(nx,1);
for i=1:nx
    x(i)=dx*i-dx;
end

% Initialisierung von Vektoren
rhs_uds=zeros(nx,1);
rhs_cds=zeros(nx,1);

% Zeitschleife
for nt=1:maxnt
    
    for i=2:nxm1
        % Konvektiver Term mit Upwind Differenzen Schema
        konv_uds=rho*u0* (phi_uds(i) - phi_uds(i-1)) / (x(i) - x(i-1));
        
        % Konvektiver Term mit zentralem Differenzen Schema
        konv_cds=rho*u0* (phi_cds(i+1) - phi_cds(i-1)) / (x(i+1) - x(i-1));
        
        % Diffusiver Term mit CDS
        diff_uds=gamma* ( (phi_uds(i+1) - phi_uds(i)) / (x(i+1) - x(i)) - (phi_uds(i) - phi_uds(i-1)) / (x(i) - x(i-1)) ) / 0.5 / (x(i+1) - x(i-1));
        diff_cds=gamma* ( (phi_cds(i+1) - phi_cds(i)) / (x(i+1) - x(i)) - (phi_cds(i) - phi_cds(i-1)) / (x(i) - x(i-1)) ) / 0.5 / (x(i+1) - x(i-1));

        rhs_uds(i)=diff_uds-konv_uds;
        rhs_cds(i)=diff_cds-konv_cds;
    end
    
    % Kalkuliere neue Werte für phi in der Zeit mit einem
    % expliziten Euler Schema
    for i=2:nxm1
        phi_uds(i)= phi_uds(i) + dt * rhs_uds(i);
        phi_cds(i)= phi_cds(i) + dt * rhs_cds(i);
    end

    % Simulationszeit
    zeit=nt*dt;

    if conv_uds == 0 && all(abs(temp_uds - phi_uds) < tol)
        conv_uds = 1;
        conv_zeit_uds = zeit;
    else
        temp_uds = phi_uds;
    end
    if conv_cds == 0 && all(abs(temp_cds - phi_cds) < tol)
        conv_cds = 1;
        conv_zeit_cds = zeit;
    else
        temp_cds = phi_cds;
    end

    % Ergebnisse werden dargestellt
    % Nach jeder n-ten Iteration soll das Ergebnis geplottet werden
    % (werden die Ergebnisse nicht in jeder Iteration dargestellt, läuft
    % die Simulation schneller)
    if rem(nt-1,n)==0
        figure(1)
        plot(xel,el,'-b',x,phi_uds,'-ro',x,phi_cds,'-kx');
        legend('exakt','uds','cds','Location','Eastoutside');
        axis([0 1 -0.2 1])
        text(0.1,0.9,['Peclet-Zahl= ',num2str(peclet),' (\rho = ',num2str(rho),', u = ',num2str(u0),', \Gamma = ',num2str(gamma),')']);
        text(0.1,0.85,['DCFL-Zahl= ',num2str(dcfl),'(CFL = ',num2str(cfl),', D = ',num2str(d),', T = ',num2str(t),')'])
        text(0.1,0.8,['Zeit= ',num2str(zeit)])
        text(0.1,0.75,['Zeitschritt= ',num2str(dt)])
        text(0.1,0.7,['Knotenpunkte= ',num2str(nx)])
        if conv_uds == 1
            text(0.1,0.6,['UDS Convergiert @ ',num2str(conv_zeit_uds),' s'])
        end
        if conv_cds == 1
            text(0.1,0.55,['CDS Convergiert @ ',num2str(conv_zeit_cds),' s'])
        end
        drawnow
        
        phi_uds_t(floor(nt/n)+1,:) = phi_uds';
        phi_cds_t(floor(nt/n)+1,:) = phi_cds';
    end
   
end
%Ende Zeitschleife

if conv_uds == 0
    text(0.1,0.6,['bis ',num2str(maxnt*dt),'s UDS Convergiert nicht'])
end
if conv_cds == 0
    text(0.1,0.55,['bis ',num2str(maxnt*dt),'s CDS Convergiert nicht'])
end

[X, Y] = meshgrid(0:dx:1,0:n*dt:dt*maxnt);
figure(2);
clf;
s1 = surface(X,Y,phi_uds_t,X);
title('\Phi of UDS wrt. time');
xlabel('x');
ylabel('t');
zlabel('\Phi');
colormap(parula);
set(s1, 'EdgeColor', 'none');
colorbar;
grid on;
view(3);

figure(3);
clf;
s2 = surface(X,Y,phi_cds_t,X);
title('\Phi of CDS wrt. time');
xlabel('x');
ylabel('t');
zlabel('\Phi');
colormap(parula);
set(s2, 'EdgeColor', 'none');
colorbar;
grid on;
view(3);

end
% Ende Funktion