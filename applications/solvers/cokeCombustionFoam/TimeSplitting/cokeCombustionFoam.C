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

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "cokeEqn.H"

        #include "rhoEqn.H"
   
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        { 
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
               #include "pEqn.H"
            }

            #include "YEqn.H"

            #include "EEqn.H"
        
        }
                
        rho = thermo.rho();
        rhoByEps=rho*rEps;

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;


    
    const volScalarField& T=thermo.T();
    scalar maxT=gMax(T);

    const volScalarField& YO2= Y[O2Index];
    scalar minO2=gMin(YO2);

    scalar maxUx=gMax(U.component(0)->field());

    volScalarField tempCoke("tempCoke",coke);
    forAll(solid,celli)
    {
        if(solid[celli]<1.0) // in the solid region
        {
            tempCoke[celli]=great;
        }
    }
    scalar minCoke=gMin(tempCoke);

    scalar endTimeSeconds=runTime.endTime().value();
    Info<<"endTime: "<<endTimeSeconds<<endl;
    
    Info<<"Final Time step: "<<runTime.deltaT().value()<<endl;
    Info<<"Max T in solid region: "<<maxT<<endl;
    Info<<"Min O2 in fluid region: "<<minO2<<endl;
    Info<<"Max Ux in fluid region: "<<maxUx<<endl;

    Info<<"Min coke in solid region: "<<minCoke<<endl;
    scalar maxBurningRate=(1-minCoke/0.8);
    Info<<"Max coke burning rate: "<<maxBurningRate*100<<"%"<<endl;
    Info<<"Max coke burning rate in one second: "<<1.0/endTimeSeconds*maxBurningRate*100<<"%"<<endl;

    
    return 0;
}