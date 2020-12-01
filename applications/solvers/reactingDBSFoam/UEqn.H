    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        (1./eps)*(fvm::ddt(rhof, U) + fvm::div(phiByEpsf, U))
      - (1.0/eps)*fvm::laplacian(nu, U)
      
     ==
        fvOptions(rhof, U) - fvm::SuSp(nu*rK,U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhof)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
