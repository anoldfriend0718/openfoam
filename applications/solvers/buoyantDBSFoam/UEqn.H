    //update the filed related to porous medium 
    phiByEpsf = phi*repsf;
    mu=thermo.mu();

    // Solve the Momentum equation
    fvVectorMatrix UEqn
    (
        rEps*(fvm::ddt(rhof, U) + fvm::div(phiByEpsf, U))
      - rEps*fvm::laplacian(mu, U)
     ==
        fvOptions(rhof, U) - fvm::SuSp(mu*rK,U)
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
