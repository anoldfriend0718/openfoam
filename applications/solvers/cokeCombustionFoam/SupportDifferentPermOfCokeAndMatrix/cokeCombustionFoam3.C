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

#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "fvcVolumeIntegrate.H"
#include "rhoReactionThermo.H"
#include "constSolidThermo.H"
#include "cokeCombustion.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "scalar.H"
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

    Info<<"Reading the porosity and coke fraction"<<nl<<endl;
    volScalarField eps
    (
        IOobject
        (
            "eps",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField coke
    (
        IOobject
        (
            "coke",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField rock
    (
        IOobject
        (
            "rock",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1-eps-coke
    );

    volScalarField rEps
    {
        "rEps",
        1.0/(eps+SMALL)
    };

    surfaceScalarField rEpsf
    {
        "rEpsf",
        fvc::interpolate(rEps)
    };

    Info<< "calculating face flux field phi in porous medium"<< nl << endl;
    surfaceScalarField phiByEpsf
    {
        "phiByEpsf",
        phi*rEpsf
    };

    volScalarField rhoByEps("rhoByEps",rho*rEps);


    Info<< "Reading porous medium transportProperties"<< nl << endl;
    IOdictionary porousTransportProperties
    (
        IOobject
        (
            "porousTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
   

    Info<< "Createing permeability field"<< nl << endl;

    dimensionedScalar K0 //rock 
    {
        "K0",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        porousTransportProperties.lookup("K0")
    };

    dimensionedScalar rK0 //rock
    {
        "rK0",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        (1.0/K0).value()
    };

    dimensionedScalar K1 //coke
    {
        "K0",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        porousTransportProperties.lookup("K1")
    };

    dimensionedScalar rK1 //coke
    {
        "rK1",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        (1.0/K1).value()
    };


    volScalarField rK
    (
        "rK",
        rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL)
    );

    //Definition of coke region indicator
    volScalarField cokeRegion
    (
        IOobject
        (
            "cokeRegion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar ("cokeRegion", dimensionSet(0,0,0,0,0,0,0), 0) 
    );

    //label coke regions 
    forAll(coke,celli)
    {
        if(coke[celli]>SMALL)
        {
            cokeRegion[celli]=1.0;
        }
        else
        {
            cokeRegion[celli]=0.0;
        }
    }

    forAll(rK,celli)
    {
        if(cokeRegion[celli]>1.0-SMALL) //==1.0
        {
            scalar epsi=eps[celli];
            rK[celli]=(1.0-epsi)*(1.0-epsi)/max(epsi*epsi*epsi,SMALL)*rK1.value();
        }
    }

    Info<<"Creating the porous medium drag force field"<< nl <<endl;
    volScalarField drag
    {
        "drag",
        fvc::average(mu*rK)
    };

    // Info<<"Label the cells that have porous medium"<<endl;
    // Definition of Solid Indicator
    volScalarField solid
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
        dimensionedScalar ("solid", dimensionSet(0,0,0,0,0,0,0), 0) 
    );


    //label solid regions 
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

    // volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
    // volScalarField::Boundary& solidBf=solid.boundaryFieldRef();
    // forAll(epsBf,patchi)
    // {
    //     forAll(epsBf[patchi],facei)
    //     {
    //         if(epsBf[patchi][facei]>0.99)
    //         {
    //             solidBf[patchi][facei]=0.0;
    //         }
    //         else
    //         {
    //             solidBf[patchi][facei]=1.0;
    //         }
    //     }
    // }

    forAll(drag,celli)
    {
        if(solid[celli]<small) //==0
        {
            drag[celli]=0.0;
        }
    }
    

    #include "initReactionModel.H"

    #include "createFvOptions.H"

    Info<< "Reading mvConvection "<< nl <<endl;
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );



    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createCombustionControl.H"

     // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run() && residualCokeFrac>targetResidualCokeFraction)
    {
        #include "readCombustionControl.H"
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "solving reaction model"<<endl;
        reaction.correct();
           
        #include "rhoEqn.H"
   
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        { 
            // #include "cokeEqn.H"

            
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

            forAll(coke,celli)
            {
                //remove the residual coke,the threshold 1e-6 defined in the specialSurfaceArea class
                if(coke[celli]<=1e-6) 
                {
                    coke[celli]=0.0;
                }
            }

            Info<<"updating the porous medium and related fields"<<endl;
            eps=1-coke-rock;
            eps.max(0.0);
            rEps=1.0/(eps+SMALL);
            rEpsf=fvc::interpolate(rEps);
            phiByEpsf=phi*rEpsf;

            //reupdate the coke region indication
            forAll(coke,celli)
            {
                if(coke[celli]>SMALL)
                {
                    cokeRegion[celli]=1.0;
                }
                else
                {
                    cokeRegion[celli]=0.0;
                }
            }

            rK=rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL);
            forAll(cokeRegion,celli)
            {
                if(cokeRegion[celli]>1.0-SMALL) //==1.0
                {
                    scalar epsi=eps[celli];
                    rK[celli]=(1.0-epsi)*(1.0-epsi)/max(epsi*epsi*epsi,SMALL)*rK1.value();
                }

            }
      
            drag=fvc::average(mu*rK);

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
            // volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
            // volScalarField::Boundary& solidBf=solid.boundaryFieldRef();
            // forAll(epsBf,patchi)
            // {
            //     forAll(epsBf[patchi],facei)
            //     {
            //         if(epsBf[patchi][facei]>0.99)
            //         {
            //             solidBf[patchi][facei]=0.0;
            //         }
            //         else
            //         {
            //             solidBf[patchi][facei]=1.0;
            //         }
            //     }
            // }

            forAll(drag,celli)
            {
                if(solid[celli]<small) //==0
                {
                    drag[celli]=0.0;
                }
            }


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

        residualCokeFrac=(fvc::domainIntegrate(coke)/gSum(mesh.V())).value();
        Info<<"Residual coke fraction: "<< residualCokeFrac<<endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;

    Info<<"Target residual coke fraction: "<<targetResidualCokeFraction<<endl;
    Info<<"Residual coke fraction: "<<residualCokeFrac<<endl;
    if(residualCokeFrac<=targetResidualCokeFraction)
    {
        Info<<"Residual coke fraction reached the target!"<<endl;
    }
    else
    {
        Info<<"Residual coke fraction did not reach the target!"<<endl;
    }
    
    const volScalarField& T=thermo.T();
    scalar maxT=gMax(T);
    Info<<"Max T: "<<maxT<<endl;

    const volScalarField& YO2= Y[O2Index];
    scalar minO2=gMin(YO2);
    Info<<"Min O2 in fluid region: "<<minO2<<endl;

    scalar maxUx=gMax(U.component(0)->field());
    Info<<"Max Ux in fluid region: "<<maxUx<<endl;

    return 0;
}