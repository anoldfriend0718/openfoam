/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    cokeCombustionFoam

Description
    Solver for combustion with chemical reactions using a density based
    thermodynamics package with enhanced buoyancy treatment.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "constSolidThermo.H"
#include "cokeCombustion.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"



int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

     // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    volScalarField rhoCoke
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoCoke", dimDensity, cokeThermo.density) 
    );


    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "solving reaction model"<<endl;
        reaction.correct();

        //solving coke volume evolution equation
        fvScalarMatrix cokeEqn
        (
            fvm::ddt(rhoCoke,coke)
          ==
              reaction.Rs(coke)
            + fvOptions(rhoCoke,coke)
        );
        cokeEqn.relax();
        fvOptions.constrain(cokeEqn);

        cokeEqn.solve();
        fvOptions.correct(coke);

        Info<<"updating the porous medium and related fields"<<endl;
        eps=1-coke-rock;
        rEps=1.0/(eps+SMALL);
        rEpsf=fvc::interpolate(rEps);
        phiByEpsf=phi*rEpsf;
        rhoByEps=rho*rEps;

        forAll(eps,celli)
        {
            if(eps[celli]>0.99)
            {
                solid[celli]=0.0;
            }
            else
            {
                solid[celli]=1.0;
            }
        }
        volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
        volScalarField::Boundary& solidBf=solid.boundaryFieldRef();
        forAll(epsBf,patchi)
        {
            forAll(epsBf[patchi],facei)
            {
                if(epsBf[patchi][facei]>0.99)
                {
                    solidBf[patchi][facei]=0.0;
                }
                else
                {
                    solidBf[patchi][facei]=1.0;
                }
            }
        }

        rK=rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL);
        drag=fvc::average(mu*rK);
        forAll(drag,celli)
        {
            if(solid[celli]<small) //==0
            {
                drag[celli]=0.0;
            }
        }





        










        // #include "cokeEqn.H"
        // #include "rhoEqn.H"

        // // --- Pressure-velocity PIMPLE corrector loop
        // while (pimple.loop())
        // { 
        //     #include "UEqn.H"

        //     #include "YEqn.H"

        //     #include "EEqn.H"

        //     // --- Pressure corrector loop
        //     while (pimple.correct())
        //     {
        //         #include "pEqn.H"
        //     }
        // }

        // rhof = thermo.rho();
        
        // #include "updateVariables.H"

        // runTime.write();
        
        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        //     << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}