
## Residual Control Hints 

- Initially, the fluid flow and heat flow is developing, leading to the difficult convergence. We can disable the reaction to calculate the steady-state fluid and pressure, and then map field to the actual case to improve the convergence at the initial stage.

```
 mapFields ../calciteCrystalWithoutReactions -consistent -sourceTime 1
```

 - At the initial stage, the energy and species equation cannot reach the thresholds of 1e-5 within 100 iterations. After the fluid flow and heat flow is quasi-steady, the convergence at each time step is quit rapid with 2-4 "pimple outer corrector" loop to get the all the residuals down to 1e-6 or even 1e-7. To consider the different convergence conditions of the initial stage and quasi-steady conditions, the  nOuterCorrectors to a very relative high number (~ 50 to 1000) and control your pimple loop with the residual control. In this case, I set the nOuterCorrectors as 100.

 ```
 PIMPLE
{
    momentumPredictor yes;
    nNonOrthogonalCorrectors 0;
    nCorrectors       2;
    nOuterCorrectors  100;

    outerCorrectorResidualControl
    {
        e     
        {
            tolerance  5e-7;
            relTol     0;
        }

        Y    
        {
            tolerance  1e-7;
            relTol     0;
        }

        U      
        {
            tolerance  1e-8;
            relTol     0;
        }

        p_rgh      
        {
            tolerance  1e-8;
            relTol     0;
        }
        
    }
 ```

- When the solid dissolution appear, the convergence of momentum and pressure equation  turn oscillation. Considering that, we need to keep a slightly tight threshold for them so that we can get smooth fluid and pressure field during the solid dissolution. But the optimum threshold was given by trial-and-error at the beginning to get the detailed understanding of the numerical behavior of the specified physical problem. One small tip: we can adjust the threshold during running and OpenFoam can take them effect immediately, since the control dict and fvSolution dict are set with MUST_READ_IF_MODIFIED read options 

- The tight residual control obviously slow down the convergence and increase the computational cost. We can implement sensitivity analysis on the dynamics of  our interested physical variables, such the solid fraction, mean/max temperature. A small tip: we can use the built-in or self-developed function object to postprocess these interested physical field, such as cellMin, cellMax.  If we want to develop our own "fvMeshFunctionObject", we can take "mag" for reference 
