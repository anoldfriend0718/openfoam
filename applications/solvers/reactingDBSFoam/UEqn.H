    //update the filed related to porous medium 
    rK = rK0*(1.0-eps)*(1.0-eps)/(eps*eps*eps);
    epsf=fvc::interpolate(eps);
    phiByEpsf = phi/epsf;
    mu=thermo.mu();

    // Solve the Momentum equation
    fvVectorMatrix UEqn
    (
        (1./eps)*(fvm::ddt(rhof, U) + fvm::div(phiByEpsf, U))
      - (1.0/eps)*fvm::laplacian(mu, U)
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
