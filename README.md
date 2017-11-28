# Modellbildung-und-Simulation
Übungen zu Modellbildung und Simulation WS17/18 in KIT

> Update 2017/11/28
> 1. Lokale Fehler sind nun mit Differenzequozienten zu rechnen.
> 2. Bei Optimierung von `h` wird zuerst die Bedingung `|x(n+1)-x(n)| > delta` beurteilt.
> 3. Euler-Implizit-Verfahren und Runge-Kutta-Verfahren mit optimierten `h` sind hinzufuegt.
> 4. Um Rechnungszeit zu sparen, Regelmäßige Matrix ist beim Bestimmen globaler Fehler benutzt.

## Uebung 6: Numerische Integration gewöhnlicher Differentialgleichungen
Integrationsverfahren:
  * Euler-explizit  
    `x(n+1) = xn + f(xn, tn)·h`
  * Euler-implizit  
    `x(n+1) = xn + f(xn+1, tn+1)·h`
  * Heun-Verfahren (2-Ordnung)  
    `k1(n) =h·f(xn, tn)`  
    `k2(n) =h·f(xn + k1(n) , tn + h)`  
    `x(n+1) =xn + 1/2 ·[k1(n) + k2(n)]`
  * Runge-Kutta-Verfahren (4-Ordnung)  
    `k1(n) =h·f(xn, tn)`  
    `k2(n) =h·f(xn + k1(n)/2 , tn + h/2)`  
    `k3(n) =h·f(xn + k2(n)/2 , tn + h/2)`  
    `k4(n) =h·f(xn + k3(n) , tn + h)`  
    `x(n+1) =xn + 1/6 ·[k1(n) + 2·k2(n) + 2·k3(n) + k4(n)]`

[MathPendel_Bearbeitungsfile.m](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/MathPendel_Bearbeitungsfile.m)
 * Systemmatrix erstellen
 * Anfangsbedingungen einsetzen
 * Lösung der Gleichung mit verschiedenen numerischen Verfahren berechnen
 * Lakale und Globale Fehler berechnen
 * Zeitschritte optimieren und vergleichen
 > <center>Numerische Verfahren mit fixen h</center>
 > ![image](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Verfahren_m_fixen_h.jpg)
 > <center>Numerische Verfahren mit optimierten h</center>
 > ![image](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Verfahren_m_optim_h.jpg)

[Schwingerkette_Bearbeitungsfile.m](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Schwingerkette_Bearbeitungsfile.m)
 * Systemmatrix erstellen
 * Anfangsbedigungen einsetzen
 * Lösung der Gleichung mit Heun-Verfahren berechnen
 * Eigenfrequenzen mit steigenden Federsteifigkeiten wieder berechnen
 * Ergebnis plotten  
 > [Schwingerkette_Model.slx](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Schwingerkette_Model.slx) stellt das o.g. Schwingersystem in `MATLAB-Simulink` auf
