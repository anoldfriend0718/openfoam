/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // #include "alphaControls.H"
            #pragma region alphaControls.H
            const dictionary& alphaControls = mesh.solverDict(alpha1.name());

            label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

            label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

            bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

            // Apply the compression correction from the previous iteration
            // Improves efficiency for steady-simulations but can only be applied
            // once the alpha field is reasonably steady, i.e. fully developed
            bool alphaApplyPrevCorr
            (
                alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
            );

            // Isotropic compression coefficient
            scalar icAlpha
            (
                alphaControls.lookupOrDefault<scalar>("icAlpha", 0)
            );

            // Shear compression coefficient
            scalar scAlpha
            (
                alphaControls.lookupOrDefault<scalar>("scAlpha", 0)
            );

            #pragma endregion alphaControls.H

            // #include "alphaEqnSubCycle.H"
            #pragma region alphaEqnSubCycle.H
            if (nAlphaSubCycles > 1)
            {
                dimensionedScalar totalDeltaT = runTime.deltaT();
                surfaceScalarField rhoPhiSum
                (
                    IOobject
                    (
                        "rhoPhiSum",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(rhoPhi.dimensions(), 0)
                );

                tmp<volScalarField> trSubDeltaT;

                if (LTS)
                {
                    trSubDeltaT =
                        fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
                }

                for
                (
                    subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                    !(++alphaSubCycle).end();
                )
                {
                    #include "alphaEqn.H"
                    rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
                }

                rhoPhi = rhoPhiSum;
            }
            else
            {
                // #include "alphaEqn.H"
                #pragma region alphaEqn.H
                {
                    word alphaScheme("div(phi,alpha)");
                    word alpharScheme("div(phirb,alpha)");

                    // Set the off-centering coefficient according to ddt scheme
                    scalar ocCoeff = 0;
                    {
                        tmp<fv::ddtScheme<scalar>> tddtAlpha
                        (
                            fv::ddtScheme<scalar>::New
                            (
                                mesh,
                                mesh.ddtScheme("ddt(alpha)")
                            )
                        );
                        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

                        if
                        (
                            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
                        || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
                        )
                        {
                            ocCoeff = 0;
                        }
                        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
                        {
                            if (nAlphaSubCycles > 1)
                            {
                                FatalErrorInFunction
                                    << "Sub-cycling is not supported "
                                    "with the CrankNicolson ddt scheme"
                                    << exit(FatalError);
                            }

                            if
                            (
                                alphaRestart
                            || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
                            )
                            {
                                ocCoeff =
                                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                                .ocCoeff();
                            }
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "Only Euler and CrankNicolson ddt schemes are supported"
                                << exit(FatalError);
                        }
                    }

                    // Set the time blending factor, 1 for Euler
                    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

                    // Standard face-flux compression coefficient
                    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

                    // Add the optional isotropic compression contribution
                    if (icAlpha > 0)
                    {
                        phic *= (1.0 - icAlpha);
                        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
                    }

                    // Add the optional shear compression contribution
                    if (scAlpha > 0)
                    {
                        phic +=
                            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
                    }


                    surfaceScalarField::Boundary& phicBf =
                        phic.boundaryFieldRef();

                    // Do not compress interface at non-coupled boundary faces
                    // (inlets, outlets etc.)
                    forAll(phic.boundaryField(), patchi)
                    {
                        fvsPatchScalarField& phicp = phicBf[patchi];

                        if (!phicp.coupled())
                        {
                            phicp == 0;
                        }
                    }

                    tmp<surfaceScalarField> phiCN(phi);

                    // Calculate the Crank-Nicolson off-centred volumetric flux
                    if (ocCoeff > 0)
                    {
                        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
                    }

                    if (MULESCorr)
                    {
                        // #include "alphaSuSp.H"
                        zeroField Su;
                        zeroField Sp;
                        zeroField divU;


                        fvScalarMatrix alpha1Eqn
                        (
                                (
                                    LTS
                                ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
                                : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
                                )
                            + fv::gaussConvectionScheme<scalar>
                                (
                                    mesh,
                                    phiCN,
                                    upwind<scalar>(mesh, phiCN)
                                ).fvmDiv(phiCN, alpha1)
                            // - fvm::Sp(fvc::ddt(dimensionedScalar(dimless, 1), mesh)
                            //           + fvc::div(phiCN), alpha1)
                            ==
                                Su + fvm::Sp(Sp + divU, alpha1)// sp should include the effect of (1/rho1-alpha1/rho1+alpha1/rho2)
                        );

                        alpha1Eqn.solve();

                        Info<< "Phase-1 volume fraction = "
                            << alpha1.weightedAverage(mesh.Vsc()).value()
                            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
                            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
                            << endl;

                        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
                        alphaPhi10 = talphaPhi1UD();

                        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
                        {
                            Info<< "Applying the previous iteration compression flux" << endl;
                            MULES::correct
                            (
                                geometricOneField(),
                                alpha1,
                                alphaPhi10,
                                talphaPhi1Corr0.ref(),
                                oneField(),
                                zeroField()
                            );

                            alphaPhi10 += talphaPhi1Corr0();
                        }

                        // Cache the upwind-flux
                        talphaPhi1Corr0 = talphaPhi1UD;

                        alpha2 = 1.0 - alpha1;

                        mixture.correct();
                    }

                    //in order to circumvent the issue of the non-linearity of the fluxes
                    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++) //make nAlphaCorr MUKES integrators 
                    {
                        // #include "alphaSuSp.H"
                        zeroField Su;
                        zeroField Sp;
                        zeroField divU;


                        surfaceScalarField phir(phic*mixture.nHatf()); //nHatf and thus phir change with acorr 

                        tmp<surfaceScalarField> talphaPhi1Un
                        (
                            fvc::flux
                            (
                                phiCN(),
                                cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime(), //alpha1 change with aCorr
                                alphaScheme
                            )
                        + fvc::flux
                            (
                            -fvc::flux(-phir, alpha2, alpharScheme),
                                alpha1,
                                alpharScheme
                            )
                        );

                        if (MULESCorr)
                        {
                            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);//alphaPhi10, the alpha flux at the previous acorr step
                            volScalarField alpha10("alpha10", alpha1);

                            MULES::correct
                            (
                                geometricOneField(),
                                alpha1,
                                talphaPhi1Un(),
                                talphaPhi1Corr.ref(),
                                Sp,
                                (-Sp*alpha1)(), //Su for the MULES::correct
                                oneField(),
                                zeroField()
                            );

                            // Under-relax the correction for all but the 1st corrector
                            if (aCorr == 0)
                            {
                                alphaPhi10 += talphaPhi1Corr();
                            }
                            else
                            {
                                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                                alphaPhi10 += 0.5*talphaPhi1Corr();
                            }
                        }
                        else
                        {
                            alphaPhi10 = talphaPhi1Un;

                            MULES::explicitSolve
                            (
                                geometricOneField(),
                                alpha1,
                                phiCN,
                                alphaPhi10,
                                Sp,
                                (Su + divU*min(alpha1(), scalar(1)))(),
                                oneField(),
                                zeroField()
                            );
                        }

                        alpha2 = 1.0 - alpha1;

                        mixture.correct();
                    }

                    if (alphaApplyPrevCorr && MULESCorr)
                    {
                        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
                        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
                    }
                    else
                    {
                        talphaPhi1Corr0.clear();
                    }

                    // #include "rhofs.H"
                    const dimensionedScalar& rho1f(rho1);
                    const dimensionedScalar& rho2f(rho2);

                    if
                    (
                        word(mesh.ddtScheme("ddt(rho,U)"))
                    == fv::EulerDdtScheme<vector>::typeName
                    || word(mesh.ddtScheme("ddt(rho,U)"))
                    == fv::localEulerDdtScheme<vector>::typeName
                    )
                    {
                        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
                    }
                    else
                    {
                        if (ocCoeff > 0)
                        {
                            // Calculate the end-of-time-step alpha flux
                            alphaPhi10 =
                                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
                        }

                        // Calculate the end-of-time-step mass flux
                        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
                    }

                    Info<< "Phase-1 volume fraction = "
                        << alpha1.weightedAverage(mesh.Vsc()).value()
                        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
                        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
                        << endl;
                }

                
                #pragma endregion alphaEqn.H
            }

            rho == alpha1*rho1 + alpha2*rho2;

            #pragma endregion alphaEqnSubCycle.H

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
