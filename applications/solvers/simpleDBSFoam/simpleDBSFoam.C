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
#include "simpleControl.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolate.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

     while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Momentum predictor

        tmp<fvVectorMatrix> tUEqn
        (
            (1./eps)*fvm::ddt(U)
          + (1./eps)*fvm::div(phiByEpsf, U)
          - (1./eps)*fvm::laplacian(nu, U)
          + fvm::SuSp(drag,U)
        );

        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();

        if (simple.momentumPredictor())
        {
            solve(UEqn == (-fvc::grad(p)));
        }

        // --- SIMPLE loop
    
        volScalarField rAU(1.0/UEqn.A());
        //add the phi caused by gravity acceleration 
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, p);

        tmp<volScalarField> rAtU(rAU);

        if (simple.consistent())
        {
            rAtU = 1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
            HbyA -= (rAU - rAtU())*fvc::grad(p);
        }

        tUEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU);

        // Non-orthogonal pressure corrector loop
        while (simple.correctNonOrthogonal())
        {
            // Pressure corrector
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        #include "continuityErrs.H"

        // Explicitly relax pressure for momentum corrector
        p.relax();

        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
        phiByEpsf = phi*repsf;


        runTime.write();

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
