// Solve the Momentum equation
Info << "solve UEqn" << endl;
MRF.correctBoundaryVelocity(U);
surfaceScalarField epsf ("epsf", fvc::interpolate(eps));
surfaceScalarField phibyeps (phi/epsf);
volScalarField rKesp ("rKesp", (1.0-eps)*(1.0-eps)/eps/eps/eps/K);   //Kozeny-Carman
/*
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
*/
tmp<fvVectorMatrix> tUEqn
(
    (1.0/eps)*fvm::ddt(U) 
  + (1.0/eps)*fvm::div(phibyeps, U)
  + MRF.DDt(U)
  - (1.0/eps)*fvm::laplacian(nu, U)
  + fvm::Sp(nu*rKesp,U)
 ==
    fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
