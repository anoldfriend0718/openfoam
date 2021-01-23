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

        //// #include "cokeEqn.H"
        Info<< "solving reaction model"<<endl;
        reaction.correct();
        cokeRectionRate=reaction.Rs(coke) & coke;
        Qdot=reaction.Qdot().ref();

        //solving coke volume evolution equation
        fvScalarMatrix cokeEqn
        (
            fvm::ddt(rhoCoke,coke)
          ==
            reaction.Rs(coke) //note signs
          + fvOptions(rhoCoke,coke)
        );
        cokeEqn.relax();
        fvOptions.constrain(cokeEqn);

        cokeEqn.solve();
        fvOptions.correct(coke);
        coke=max(min(coke,1.0),0.0);

        Info<<"updating the porous medium and related fields"<<endl;
        eps=1-coke-rock;
        rEps=1.0/(eps+SMALL);
        rEpsf=fvc::interpolate(rEps);
        phiByEpsf=phi*rEpsf;

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

            //     #include "EEqn.H"
            {
                volScalarField& he = thermo.he();

                cokeRhoCpByCpvf.storeOldTimes();
                cokeRhoCpByCpvf=cokeThermo.density*cokeThermo.Cp/thermo.Cpv();
                
                rockRhoCpByCpvf.storeOldTimes();
                rockRhoCpByCpvf=rockThermo.density*rockThermo.Cp/thermo.Cpv();

                tmp<volScalarField> talphaEff = eps * thermo.alphahe()
                                              + coke* cokeThermo.kappa/thermo.Cpv()
                                              + rock* rockThermo.kappa/thermo.Cpv();
                const volScalarField& alphaEff = talphaEff();    
                
                //only suitable for the const cpv    
                fvScalarMatrix EEqn
                (
                      fvm::ddt(eps, rho,            he) 
                    + fvm::ddt(coke,cokeRhoCpByCpvf,he)
                    + fvm::ddt(rock,rockRhoCpByCpvf,he)
                    + mvConvection->fvmDiv(phi, he)
                    + fvc::ddt(rho, K) + fvc::div(phi, K)
                    + (
                          he.name() == "e"
                        ? fvc::div
                            (
                                fvc::absolute(phi/fvc::interpolate(rho), U),
                                p,
                                "div(phiv,p)"
                            )
                        : -dpdt
                      )
                    - fvm::laplacian(alphaEff, he)
                  ==
                      rho*(U&g)
                    + fvm::Su(reaction.Qdot(),he)
                    + fvOptions(rho, he)
                );   

                EEqn.relax();

                fvOptions.constrain(EEqn);

                EEqn.solve();

                fvOptions.correct(he);

                thermo.correct();

                const volScalarField& T=thermo.T();

                Info<< "min/max(T) = "
                    << min(T).value() << ", " << max(T).value() << endl;    
            }
        }
                
        rho = thermo.rho();
        rhoByEps=rho*rEps;

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}