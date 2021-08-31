// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
{
      rEps*(fvm::ddt(rho, U) + fvm::div(phiByEpsf, U))
    - rEps*fvm::laplacian(mu, U) //mu is the reference of thermo.mu()
  ==
    - fvm::SuSp(drag,U)
    + fvOptions(rho, U) 
};
fvVectorMatrix& UEqn = tUEqn.ref();

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
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
