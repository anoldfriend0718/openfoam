{
    // Info<<"solving pressure equation"<<endl;

    rho = thermo.rho();
    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
             fvc::flux(rho*HbyA)
            //rhofByEps old time is called in ddtCorr,
            // so we need to update the rhofByEps 
            // when one time step completes or a new time step starts
            + rhorAUf*fvc::ddtCorr(rhoByEps, U, phiByEpsf)  
        )
    );

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);
    // Info<<"closedVolume:"<<closedVolume<<endl;

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf);
    // fvScalarMatrix p_rghDDtEqn
    // (
    //     eps*fvc::ddt(rho) 
    //     + eps*psi*correction(fvm::ddt(p_rgh))
    //     + fvc::div(phiHbyA)
    // ==
    //     fvm::Su(RRg,p_rgh)
    //     + fvOptions(psi, p_rgh, rho.name())
    // );


    //update the RRg source term due to the rho change after correcting thermo model 
    tmp<volScalarField> tRRg=volScalarField::New
    (
        "zeroRRg",
        mesh,
        dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
    );

    if(involveGasProduction)
    {
        const tmp<volScalarField>& deltarRho=1./rho-1.0/rhoCoke;
        tRRg=(-deltarRho*rho*reaction.Rs(coke))&coke; //note signs
    }
    const volScalarField& RRg=tRRg.ref();

    fvScalarMatrix p_rghDDtEqn
    (
          eps*fvc::ddt(rho) 
        + eps*psi*correction(fvm::ddt(p_rgh))
        + fvc::div(phiHbyA)
      ==
        fvm::Su(RRg,p_rgh)
        + fvOptions(psi, p_rgh, rho.name())
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

    p = p_rgh + rho*gh;

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (!thermo.incompressible() && closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh;
    }
    // Thermodynamic density update
    thermo.correctRho(psi*p - psip0);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

}