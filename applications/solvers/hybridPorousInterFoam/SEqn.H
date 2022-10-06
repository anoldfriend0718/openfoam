        // volScalarField mdotsi ("mdotsi", rhof*gama*alpha_si*eps*eps_si*(C-Ceq)); 
        Info << "solve SEqn" << endl;
        //volScalarField mdots ("mdots", mdotc*gama_si/rho_si); 
        //mdots = mdotc*gama_si/rho_si;       
        Info << "mdots" << " Min = " << gMin(mdots) << " Max = " << gMax(mdots) << endl;
        fvScalarMatrix SEqn
        (
            fvm::ddt(eps_si)   
          // - fvm::Sp(divU, C)    // add 0803
            ==
          //- fvm::Sp(mdotf*gama_si/rhos,eps_si)
          // - fvm::Su(mdotf*gama_si/rhos,eps_si)
            fvm::Su(mdots,eps_si)
        );

        SEqn.relax();
        //fvOptions.constrain(SEqn);
        SEqn.solve();
        //fvOptions.correct(eps_si);
        Info << "SEqn" << " Min(eps_si) = " << gMin(eps_si) << " Max(eps_si) = " << gMax(eps_si) << endl;
        eps_si = max(min(eps_si,0.999),0.0);
        eps = 1.0-eps_s0-eps_si;
        Info << "limit: SEqn" << " Min(eps_si) = " << gMin(eps_si) << " Max(eps_si) = " << gMax(eps_si) << endl;
        Info << "eps" << " Min(eps) = " << gMin(eps) << " Max(eps) = " << gMax(eps) << endl;
        eps = max(min(eps,1.0),0.001);
        Info << "limit: eps" << " Min(eps) = " << gMin(eps) << " Max(eps) = " << gMax(eps) << endl;
