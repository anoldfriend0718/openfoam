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
        //rhofByEps old time is called in ddtCorr, so we need to update the rhofByEps when one time step completes or a new time step starts
      + rhorAUf*fvc::ddtCorr(rhofByEps, U, phiByEpsf)  
    )
);

bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);
// Info<<"closedVolume:"<<closedVolume<<endl;

phiHbyA += phig;
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

// For closed-volume cases adjust the pressure level
// to obey overall mass continuity
if (!thermo.incompressible() && closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
    p_rgh = p - rhof*gh;
}


// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"
