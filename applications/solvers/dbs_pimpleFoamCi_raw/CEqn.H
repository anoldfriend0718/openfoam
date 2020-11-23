        Info << "solve CEqn" << endl;
        volScalarField eps2D (eps*eps*Df);
        volScalarField mdotf ("mdotf", mag(fvc::grad(eps))*(4.0*eps*(1-eps))*r_si*rhof); 
        // volScalarField mdotf ("mdotf", rhof*alpha_si*eps*eps_si);  
        Info << "mdotf" << " Min = " << gMin(mdotf) << " Max = " << gMax(mdotf) << endl;  
        fvScalarMatrix CEqn
        (
            fvm::ddt(eps, C)   
          + fvm::div(phi, C)
          - fvm::laplacian(eps2D, C)
          // - fvm::Sp(divU, C)    // add 0803
             ==
          - fvm::Sp(mdotf/rhof,C)
          + fvm::Su(mdotf/rhof*Ceq,C)
          // - fvm::Su(mdotf/rhof,C)
        );

        CEqn.relax();
        fvOptions.constrain(CEqn);
        CEqn.solve();
        fvOptions.correct(C);
        Info << "CEqn" << " Min(C) = " << gMin(C) << " Max(C) = " << gMax(C) << endl;
        C = max(min(C,1.0),0.0);
        Info << "limit: CEqn" << " Min(C) = " << gMin(C) << " Max(C) = " << gMax(C) << endl;
        
