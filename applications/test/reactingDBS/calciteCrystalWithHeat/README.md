
## Residual Control Hints 

- Initially, the fluid flow and heat flow is developing, leading to the difficult convergence. At this moment, the energy and species equation cannot reach the thresholds of 1e-5 within 100 iterations. After the fluid flow and heat flow is quasi-steady, the convergence at each time step is quit rapid with 2-4 "pimple outer corrector" loop to get the all the residuals down to 1e-6. 

- To consider the different convergence conditions of the initial stage and quasi-steady conditions, the  nOuterCorrectors to a very relative high number (~ 50 to 1000) and control your pimple loop with the residual control (make sure that explizit terms converged). In this case, I set the nOuterCorrectors as 100 

- When the solid dissolution appear, the convergence of momentum and pressure equation slightly turn worse. Considering that, we need to keep a tight threshold for them so that we can get smooth fluid and pressure field during the solid dissolution 