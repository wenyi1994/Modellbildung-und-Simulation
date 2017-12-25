function [pp, up, vp] = MuS_FDM_2d(varargin)
% Modellbildung und Simulation
% Übung zum Kapitel: Zeitkontinuierliche Modelle mit verteilten Parametern
% Balazs Pritz pritz@kit.edu
% 2D Testfall: Kavität mit drei stehenden Wänden und mit einer bewegenden Wand
% Diskretisierungsschema: Finite Differenzen Methode
% Ableitungen gerechnet mit: zentralem Differenzen Schema
% Quelle: Griebel, Michael: Numerische Simulation in der Strömungsmechanik,
% Vieweg, 1995
% ----------------------------------------------------------------------
% [pp, up, vp] = MuS_FDM_2d('dt',dt,'nt',nt,'u0',u0,'nx',nx,'ny',ny,'rho',rho,'mue',mue,'gamma',gamma,'output',n_output)
% [pp, up, vp] = MuS_FDM_2d('default')
% [pp, up, vp] = MuS_FDM_2d('default','dt',dt,'nt',nt,'material','air', ...)
% ----------------------------------------------------------------------
% input:
% dt                    Zeitschritt
% maxnt                 Maximale Anzahl von Zeitschritten
%                       (Simulationszeit=timestep*maxnt)
%                       (Es ist auch möglich, dass die Simulation nicht nur
%                       durch maxnt beendet wird, sondern die Änderung der
%                       Lösung von zwei Zeitschritten berechnet wird und
%                       anhand dieser ein Abbruchkriterium definiert wird.)
% u0                    Geschwindigkeit des oberen Randes (y=1)
% nx                    Anzahl der Knotenpunte
% ny                    Anzahl der Knotenpunte
% rho                   Dichte, Anzahl oder 'water', 'air', ...
% mue                   Dynamische Viskosität, Anzahl oder 'water', 'air', ...
% material              Materialtyp, 'water', 'air', ...
%                       Der Parameter kann durch 'rho' oder 'mue' übereinander
%                       geschrieben werden.
% gamma                 Die konvektieve Terme werden mit einer Kombination
%                       von zentralen Differenzen und upwind-Differenzen
%                       gerechnet. Die Gewichtung wird durch gamma geregelt
%                       (gamma ist hier kein Diffusionskoeffizient!).
%                       gamma=1 -> reiner upwind Diff.
%                       gamma=0 -> reiner zentrale Diff.
% n_output              Anzahl der insgesamt gegebenen Grafiken
% ----------------------------------------------------------------------
% output:
% phi_uds               Lösung der PDGL mit UDS-Verfahren
% phi_cds               Lösung der PDGL mit CDS-Verfahren
% ----------------------------------------------------------------------
% Wen Yi, Karlsruhe Institut of Technology
% yi.wen@student.kit.edu
% 2017/12/25

% Parametern
if strcmp(varargin{1},'default')
    % Default Setting
    dt = 0.1;
    maxnt = 400;
    u0 = 10.0;
    nx = 21;
    ny = 21;
    rho = 1.25;
    mue = 0.00001808;
    gamma = 0.8;
    n_output = 40;
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
        case 'ny'
        ny = varargin{i+1};
        case 'rho'
        switch varargin{i+1}
            case 'air'
            rho = 1.25;
            case 'water'
            rho = 998;
            otherwise
            rho = varargin{i+1};
        end
        case 'mue'
        switch varargin{i+1}
            case 'air'
            mue = 0.00001808;
            case 'water'
            mue = 0.001895;
            otherwise
            mue = varargin{i+1};
        end
        case 'material'
        switch varargin{i+1}
            case 'air'
            rho = 1.25;
            mue = 0.00001808;
            case 'water'
            rho = 998;
            mue = 0.001895;
            otherwise
            warning(['unexpected material ',varargin{i+1}])
            return;
        end
        case 'gamma'
        gamma = varargin{i+1};
        case 'output'
        n_output = varargin{i+1};
    end
end

% Räumliche Auflösung
% Problemgebiet: x=[0,1], y=[0,1]
% Indizierung: die Wände sind bei Index i=2, i=nx+1 bzw. j=2 und j=ny+1
% der erste Knotenpunkt in der Strömung bei Index i=3, i=nx, j=3, j=ny
% Dummy-Knoten (1 bzw. nx+2) werden benutzt, um Randbedingungen
% für die Geschwindigkeit zu definieren
%
% Es wird hier ein äquidistantes Netz verwendet.
% Wenn Sie eine lokale Verfeinerung benutzen möchten, dann müssen Sie
% x(i,j) und y(i,j) definieren, und in den Gleichungen dx und dy durch
% x(i)-x(i-1), x(i+1)-x(i) oder x(i+1)-x(i-1) ersetzen
% (natürlich dy durch y()-y()).
nxm1=nx-1;
nym1=ny-1;
nxp1=nx+1;
nyp1=ny+1;
nxp2=nx+2;
nyp2=ny+2;


dx=1/(nxm1);
dy=1/(nym1);

nue=mue/rho;

% Initialisierung
u=zeros(nxp2,nyp2);
v=zeros(nxp2,nyp2);
p=zeros(nxp2,nyp2);

% Koeffizienten für die Ableitungen in der Randnähe
ew=ones(nxp1,1);
eo=ones(nxp1,1);
es=ones(nyp1,1);
en=ones(nyp1,1);
ew(2)=0;
eo(nxp1)=0;
es(2)=0;
en(nyp1)=0;

rhs=zeros(nxp2,nyp2);

% Zeitschleife
for nt=1:maxnt
    
    % Randbedingungen für die Geschwindigkeiten
    for i=2:nxp1
        u(i,1)=-u(i,2);
        u(i,nyp2)=2.*u0-u(i,nyp1);
        v(i,1)=0;
        v(i,nyp1)=0;
    end
    for j=2:nyp1
        u(1,j)=0;
        u(nxp1,j)=0;
        v(1,j)=-v(2,j);
        v(nxp2,j)=-v(nxp1,j);
    end

    % konvektive und viskose Terme
    F=zeros(nxp1,nyp1);
    G=zeros(nxp1,nyp1);
    for j=2:nyp1
        for i=2:nx
            %
            duudx=0.25*( (u(i,j)+u(i+1,j))^2 - (u(i-1,j)+u(i,j))^2 + gamma*( abs(u(i,j)+u(i+1,j))*(u(i,j)-u(i+1,j)) - abs(u(i-1,j)+u(i,j))*(u(i-1,j)-u(i,j)) ) )/dx;
            duvdy=0.25*( (v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1)) - (v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j)) + gamma*( abs(v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1)) - abs(v(i,j-1)+v(i+1,j-1))*(u(i,j-1)-u(i,j)) ) )/dy;
            d2udx2=(u(i+1,j)-2.*u(i,j)+u(i-1,j))/dx^2;
            d2udy2=(u(i,j+1)-2.*u(i,j)+u(i,j-1))/dy^2;
            %
            % Term für die zeitliche Diskretisierung
            F(i,j)=dt*(nue*(d2udx2+d2udy2) - duudx - duvdy) + u(i,j);
        end
    end
    for j=2:ny
        for i=2:nxp1
            %
            duvdx=0.25*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) - (u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j)) + gamma*( abs(u(i,j)+u(i,j+1))*(v(i,j)-v(i+1,j)) - abs(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)-v(i,j)) ) )/dx;
            dvvdy=0.25*( (v(i,j)+v(i,j+1))^2 - (v(i,j-1)+v(i,j))^2 + gamma*( abs(v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1)) - abs(v(i,j-1)+v(i,j))*(v(i,j-1)-v(i,j)) ) )/dy;
            d2vdx2=(v(i+1,j)-2.*v(i,j)+v(i-1,j))/dx^2;
            d2vdy2=(v(i,j+1)-2.*v(i,j)+v(i,j-1))/dy^2;
            %
            % Term für die zeitliche Diskretisierung
            G(i,j)=dt*(nue*(d2vdx2+d2vdy2) - duvdx - dvvdy) + v(i,j);
        end
    end

    
    % Poissongleichung für den Druck wird gelöst
    
    % Der Relaxationsfaktor soll zwischen 0 und 2 gewählt werden
    % in der Praxis wird oft 1.7 verwendet, für 1 ergibt sich
    % das Gauß-Seidel-Verfahren
    omega=1.7;
    
    % maximale Anzahl von Iterationen bei Lösung der Poissongleichung
    pitmax=100;
    
    for pit=1:pitmax
        p_old=p;

        %err=zeros(nxp1,nyp1);

        % Randwerte werden aktualisiert
        for i=2:nxp1
            p(i,1)=p(i,2);
            p(i,nyp2)=p(i,nyp1);
            G(i,1)=v(i,1);
            G(i,nyp1)=v(i,nyp1);
        end
        for j=2:nyp1
            p(1,j)=p(2,j);
            p(nxp2,j)=p(nxp1,j);
            F(1,j)=u(1,j);
            F(nxp1,j)=u(nxp1,j);
        end

        for j=2:nyp1
            for i=2:nxp1
                rhs(i,j)=((F(i,j)-F(i-1,j))/dx + (G(i,j)-G(i,j-1))/dy)/dt;
                p(i,j)=(1.-omega)*p_old(i,j) + omega*((eo(i)*p_old(i+1,j)+ew(i)*p(i-1,j))/dx^2 + (en(j)*p_old(i,j+1)+es(j)*p(i,j-1))/dy^2 - rhs(i,j))/((eo(i)+ew(i))/dx^2 + (en(j)+es(j))/dy^2);
%                err(i,j)=(eo(i)*(p(i+1,j)-p(i,j)) - ew(i)*(p(i,j)-p(i-1,j)))/dx^2 + (en(j)*(p(i,j+1)-p(i,j)) - es(j)*(p(i,j)-p(i,j-1)))/dy^2 - rhs(i,j);
            end
        end
        % Wie die Lösung der Poissongleichung konvergiert, kann mit Hilfe
        % von err bzw. errit beobachtet werden
%        errit=0;
%        for j=2:nyp1
%            for i=2:nxp1
%                errit=errit+err(i,j)^2;
%            end
%        end
%        errit=errit/(nxp1*nyp1);
%        % errit ist der aufsummierte Fehler der aktuellen Iteration
%        % Hier könnten Sie noch ein Abbruchkriterium definieren.
%        if errit < 
%        end
    end
    
    % Die neuen Werte für die Geschwindigkeiten in der Zeit
    % wird hier mit dem expliziten Euler Schema gerechnet
    for j=2:nyp1
        for i=2:nx
            
            dpdx=(p(i+1,j)-p(i,j))/dx;
            
            % Term für die zeitliche Diskretizierung
            u(i,j)=F(i,j) - dt*dpdx/rho;
            
        end
    end
    for j=2:ny
        for i=2:nxp1
            
            dpdy=(p(i,j+1)-p(i,j))/dy;
            
            % Term für die zeitliche Diskretizierung
            v(i,j)=G(i,j) - dt*dpdy/rho;
            
        end
    end
    
    % Wenn die CFL-Zahl zu groß wird, kann die Rechnung der zeitlichen
    % Ableitung instabil werden.
    % In so einem Fall müssen Sie den Zeitschritt reduzieren.
    dcfl=0;
    for j=2:ny
        for i=2:nxp1
            
            dcfl=max(2*nue*dt/(dx*dx) + max(u(i,j)/dx,v(i,j)/dy)*dt,dcfl);
            
        end
    end
     
    % Simulationszeit
    zeit=nt*dt;

    % Output
    % Lösungsmatrizen werden für Output transformiert
    % Nach der n-ten Iteration soll das Ergebnis geplottet werden
    % (werden die Ergebnisse nicht in jeder Iteration dargestellt, läuft
    % die Simulation schneller)
    % Die Druckverteilung wird als Konturplot dargestellt mit
    % Geschwindigkeitsvektoren.
    n=floor(maxnt/n_output);
    if rem(nt,n)==0
        % Randbedingungen werden für Plot aktualisiert
        for i=2:nxp1
            u(i,1)=-u(i,2);
            u(i,nyp2)=2.*u0-u(i,nxp1);
            v(i,1)=0;
            v(i,nyp1)=0;
        end
        for j=2:nyp1
            u(1,j)=0;
            u(nxp1,j)=0;
            v(1,j)=-v(2,j);
            v(nxp2,j)=-v(nxp1,j);
        end
    
        pp=zeros(ny,nx);
        up=zeros(ny,nx);
        vp=zeros(ny,nx);
        for j=1:ny
            for i=1:nx
                pp(j,i)=p(i+1,j+1);
                up(j,i)=u(i+1,j+1);
                vp(j,i)=v(i+1,j+1);
            end
        end
        figure(1)
        [X,Y] = meshgrid(0:dx:1,0:dy:1);
        contourf(X,Y,pp)
        colorbar
        hold on
        quiver(X,Y,up,vp,'w')
        text(0,1.05,['DCFL= ',num2str(dcfl)])
        text(0.3,1.05,['Zeit= ',num2str(zeit)])
        hold off
        drawnow
    end
end
%Ende Zeitschleife

end
% Ende Funktion

