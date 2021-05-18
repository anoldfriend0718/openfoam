/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    icoGravityFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids with 
    the user-specified gravity acceleration

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvmSup.H"
#include "pisoControl.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolate.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"

    Info<< "\nReading structures" << endl;
    volScalarField solid
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    IOdictionary rockProperties
    (
        IOobject
        (
            "rockProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rockPorosity
    (
        "rockPorosity",
        dimless,
        rockProperties.lookup("porosity")
    );

    dimensionedScalar rockPerm
    (
        "rockPerm",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        rockProperties.lookup("K")
    );

    volScalarField eps
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("eps", dimless, 0) 
    );

    volScalarField rK
    (
        IOobject
        (
            "rK",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rk", dimensionSet(0, -2, 0, 0, 0, 0, 0), 0) 
    );


    forAll(solid,celli)
    {
        if(solid[celli]==0) //fluid zone 
        {
            eps[celli]=1.0;
            rK[celli]=0.0;
        }
        else //rock zone 
        {
            eps[celli]=rockPorosity.value();
            rK[celli]=1.0/(rockPerm.value());
        }
    }

    forAll(mesh.boundary(), patchi) 
    {
        forAll(solid.boundaryField()[patchi],facei)
        {
            if(solid.boundaryField()[patchi][facei]==0)   //fluid zone 
            {
                eps.boundaryFieldRef()[patchi][facei]=1.0;
                rK.boundaryFieldRef()[patchi][facei]=0.0;
            }
            else //rock zone 
            {
                eps.boundaryFieldRef()[patchi][facei]=rockPorosity.value();
                rK.boundaryFieldRef()[patchi][facei]=1.0/(rockPerm.value());
            }
        }
    }

    // eps.write();
    // rK.write();

    //Drag Coefficient Calculation
    volScalarField drag ("drag", nu*rK);

    forAll(drag,celli)
    {
        if (solid[celli]==0)//fluid zone 
        {
            drag[celli]= 0;
        }
    }

    volScalarField rEps
    (
        "rEps",
        1./eps
    );

    surfaceScalarField rEpsf
    (
        "rEpsf",
        fvc::interpolate(1/eps)
    );

    surfaceScalarField phiByEpsf
    (
        "phiByEpsf",
        phi*rEpsf
    );


    const label patchIdInlet = mesh.boundary().findPatchID("inlet");
    const fvPatch& patchInlet = mesh.boundary()[patchIdInlet];
    scalarField pInlet(patchInlet.size(),0.0);

    forAll(patchInlet,facei)
    {
        pInlet[facei]=p.boundaryFieldRef()[patchIdInlet][facei];
    }
    scalar meanPInlet=gAverage(pInlet);


    const label patchIdOutlet = mesh.boundary().findPatchID("outlet");
    const fvPatch& patchOutlet = mesh.boundary()[patchIdOutlet];
    scalarField pOutlet(patchOutlet.size(),0.0);

    forAll(patchOutlet,facei)
    {
        pOutlet[facei]=p.boundaryFieldRef()[patchIdOutlet][facei];
    }
    scalar meanPOutlet=gAverage(pOutlet);


    Info<<"average pressure at the inlet: "<<meanPInlet
        <<"; average pressure at the outlet: "<<meanPOutlet
        <<endl;

    scalar L=0.02;  //need read 
    scalar gradPressure=(meanPInlet-meanPOutlet)/L;
    Info<<"gradient of pressure: "<<gradPressure<<endl;
    scalar perm=0.0;
    


    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            rEps*fvm::ddt(U)
          + rEps*fvm::div(phiByEpsf, U)
          - rEps*fvm::laplacian(nu, U)
          + fvm::SuSp(drag,U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == (-fvc::grad(p) + g ));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            //add the phi caused by gravity acceleration 
            surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phig(rAUf*(g & mesh.Sf()));

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
                // + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                + fvc::interpolate(rAU)*fvc::ddtCorr(rEps,U, phiByEpsf)
            );

            adjustPhi(phiHbyA, U, p);

            phiHbyA+=phig;

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"


            U = HbyA - rAU*fvc::grad(p)+rAU*g;
            U.correctBoundaryConditions();
            phiByEpsf = phi*rEpsf;
        }

        runTime.write();

        // gAverage(U.component(0).)
        tmp<volScalarField> tUx=U.component(0);
        const volScalarField& Ux=tUx.ref();
        const scalar meanUx=gAverage(Ux);
        Info<<"average Ux: "<<meanUx<<endl;
        perm=meanUx*nu.value()/gradPressure*1e15; //mD
        Info<<"Permeability: "<<perm<< " mD"<<endl;
        




        tmp<volVectorField> tUResidual(U-U.oldTime());
        const volVectorField& UResidual=tUResidual.ref();
        const Foam::Vector<scalar> normUResidual=gSumCmptProd(UResidual,UResidual);
        const Foam::Vector<scalar> normU=gSumCmptProd(U,U);
        Foam::Vector<scalar> relativeUResidual;
        for (direction cmpt=0; cmpt < normU.size(); cmpt++)
        {
            relativeUResidual[cmpt] = Foam::sqrt(normUResidual[cmpt]/(normU[cmpt]+SMALL));
        }
        Info<<"relative error: Ux: "<<relativeUResidual[0]<<" Uy: "<<relativeUResidual[1]<<endl;



        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
