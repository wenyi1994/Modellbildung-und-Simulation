# Modellbildung-und-Simulation `2017/12/24`
Übungen zu Modellbildung und Simulation WS17/18 in KIT

## Uebung 10: Modelle mit verteilten Parametern - Strömungssimulation mittels FDM
### 1D Testfall: Konvektions-Diffusionsgleichung
[MuS_FDM_1d.m](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung10/MuS_FDM_1d.m)
* Räumliche Ableitung für 1. Ableitung: UDS (Upwind Difference Scheme) und CDS (Central Difference Scheme).
* Räumliche Ableitung für 2. Ableitung: CDS.
* Die Lösung wird iterativ gesucht. Es wird der Zeitterm zu der oben angegebenen Gleichung addiert. Hier wird explizite Euler Schema implementiert.  
> *Output Grafik mit `'default'` Parametern* 
> ![image](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung10/1_default.jpg)

> **Update 2017/11/29**  
> Bei Optimierung von `h` geht die Laufzeit von Anfangszeit aus zurück, um vorgegebenen Fehler zu beschränken.

> **Update 2017/11/28**  
> 1. Lokale Fehler sind nun mit Differenzequozienten zu rechnen.
> 2. Bei Optimierung von `h` wird zuerst die Bedingung `|x(n+1)-x(n)| > delta` beurteilt.
> 3. Euler-Implizit-Verfahren und Runge-Kutta-Verfahren mit optimierten `h` sind hinzufügt.
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
 > *Numerische Verfahren mit fixen h*
 > ![image](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Verfahren_m_fixen_h.jpg)
 > *Numerische Verfahren mit optimierten h*
 > ![image](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Verfahren_m_optim_h.jpg)

[Schwingerkette_Bearbeitungsfile.m](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Schwingerkette_Bearbeitungsfile.m)
 * Systemmatrix erstellen
 * Anfangsbedigungen einsetzen
 * Lösung der Gleichung mit Heun-Verfahren berechnen
 * Eigenfrequenzen mit steigenden Federsteifigkeiten wieder berechnen
 * Ergebnis plotten  
 > [Schwingerkette_Model.slx](https://github.com/wenyi1994/Modellbildung-und-Simulation/blob/master/Uebung6/Schwingerkette_Model.slx) stellt das o.g. Schwingersystem in `MATLAB-Simulink` auf
