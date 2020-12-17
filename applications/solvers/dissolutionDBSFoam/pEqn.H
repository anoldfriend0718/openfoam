rhof = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhof*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rhof)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rhof*HbyA)
      + rhorAUf*fvc::ddtCorr(rhof, U, phi)
    )
  + phig
);


// Update the pressure BCs to ensure flux consistency

constrainPressure(p_rgh, rhof, U, phiHbyA, rhorAUf);

// multiple the transient term by eps to ensure the DBS continuity equation

mdotSolid=ssArea*reactionRate*rhof*(Y-Yeq)*beta;

fvScalarMatrix p_rghDDtEqn
(
    eps*fvc::ddt(rhof) + eps*psi*correction(fvm::ddt(p_rgh))
  + fvc::div(phiHbyA)
  ==
    fvm::Su(mdotSolid*(1.0-rhof/rhos),p_rgh)
  + fvOptions(psi, p_rgh, rhof.name())
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        p_rghDDtEqn
      - fvm::laplacian(rhorAUf, p_rgh)
    );

    if (thermo.incompressible())
    {
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    }

    p_rghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA + p_rghEqn.flux(); 

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
}

p = p_rgh + rhof*gh;

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"