# Modellbildung-und-Simulation
Übungen zu Modellbildung und Simulation WS17/18 in KIT

## Uebung 6: Numerische Integration gewöhnlicher Differentialgleichungen
Integrationsverfahren:
  * Euler-explizit  
    `x(n+1) = xn + f(xn, tn)·h`
  * Euler-implizit  
    `x(n+1) = xn + f(xn+1, tn+1)·h`
  * Heun-Verfahren (2-Ordnung)  
    `k1(n) =h·f(xn, tn)`  
    `k2(n) =h·f(xn + k(1)n , tn + h)`  
    `x(n+1) =xn + 1/2 ·(k1(n) + k2(n))`
