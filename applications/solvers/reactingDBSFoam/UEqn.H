
//update the filed related to porous medium 
tEps=max(eps,minEps);
rK= rK0*(1.0-eps)*(1.0-eps)/(tEps*tEps*tEps);
epsf=fvc::interpolate(eps);
phiByEpsf = phi/epsf;
mu=thermo.mu();

// Solve the Momentum equation
fvVectorMatrix UEqn
(
      (1.0/eps)*(fvm::ddt(rhof, U) + fvm::div(phiByEpsf, U))
    - (1.0/eps)*fvm::laplacian(mu, U)
    ==
    - fvm::SuSp(mu*rK,U)
    + fvOptions(rhof, U) 
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
