    volScalarField rAU("rAU", 1.0/UEqn.A()); 
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); 
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    //volScalarField mdotp ("mdotp", mdotc*(gama_si/rho_si+gama_l/rho1+gama_g/rho2));  //liuzy, define reaction source term to mas conservation equaion
    // mdotp = mdotc*(gama_si/rho_si+gama_l/rho1+gama_g/rho2);

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    adjustPhi(phiHbyA, U, p);

    // Gravity Contribution
    surfaceScalarField phiG 
    (
       rAUf*fvc::interpolate(rho)*(g & mesh.Sf())
    );

    // Interfacial Force Contributions
    surfaceScalarField phiPc
    (
        (
            mixture.surfaceTensionForce()*Solidf0 //*Solidf0 //*(1.-Solidf) //(1.-Solidf)  //Solidf0 //(1.-Solidf) //activate free fluid interface model
	 +  fvc::interpolate(PcCoeff)*fvc::snGrad(alpha1)*(1.-Solidf0)  //*(1.-Solidf0) //*Solidf //Solidf //(1.-Solidf0)   //*Solidf //activate porous media model
        )*rAUf*mesh.magSf()
    );

    phiHbyA += (phiPc + phiG);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA) - fvm::Su(mdotp, p)   //liuzy, add reaction source term to mas conservation equaion    HbyA=laplacain+m => laplacain=HbyA-m    
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        pEqn.solve();
        // p = p*Solidv;
        if (pimple.finalNonOrthogonalIter())
        {
            phi = (phiHbyA - pEqn.flux());  //*(1.0-Solidwf)

            p.relax();

            U = HbyA + rAU*fvc::reconstruct(((phiPc+phiG) - pEqn.flux())/rAUf);
            U.correctBoundaryConditions();
	    #include "correctUBc.H"
        }
    }

    #include "continuityErrs.H"

    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }



    

