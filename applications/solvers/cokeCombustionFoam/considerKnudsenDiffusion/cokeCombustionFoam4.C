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
#include "fvcVolumeIntegrate.H"
#include "rhoReactionThermo.H"
#include "constSolidThermo.H"
#include "cokeCombustion.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "constant.H"



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
    #include "createCombustionControl.H"


    volScalarField DeffO2
    (
        IOobject
        (
            "Deff_O2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    );

    volScalarField DeffCO2
    (
        IOobject
        (
            "Deff_CO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    );




    

    






    // Info<<"M O2: "<<MO2<<endl;
    // Info<<"M CO2: "<<MCO2<<endl;
    // Info<<"M N2: "<<MN2<<endl;


    // volScalarField DO2CO2
    // (
    //     IOobject
    //     (
    //         "DO2_CO2",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    // );

    // volScalarField DO2N2
    // (
    //     IOobject
    //     (
    //         "DO2N2",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    // );

    // volScalarField DO2b
    // (
    //     IOobject
    //     (
    //         "DO2_b",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    // );

    // volScalarField DO2k
    // (
    //     IOobject
    //     (
    //         "DO2_k",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimensionSet(0,2,-1,0,0,0,0),0.0)
    // );



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
            #include "cokeEqn.H"

            {
                const volScalarField& mu=thermo.mu()();
                const volScalarField& T=thermo.T();
                
                forAll(cokeRegion,i)
                {
                    if(cokeRegion[i]>1-SMALL)
                    {
                        scalar mui=mu[i];
                        scalar pi=p[i];
                        scalar Ti=T[i];
                        scalar kni=2*(mui/pi*std::sqrt(PI*Rg1000*Ti/(2*MN2)))/dp;
                        scalar alphai=128/(15*PI2)*std::atan(4*std::pow(kni,0.4));
                        scalar rKnanoi=1./((1+alphai*kni)*(1+4*kni/(1+kni)));
                        rKnano[i]=rKnanoi;
                    }
                    else
                    {
                        rKnano[i]=1.0;
                    }
                }

                drag=fvc::average(mu*rK*rKnano);

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

                forAll(drag,celli)
                {
                    if(solid[celli]<small) //==0
                    {
                        drag[celli]=0.0;
                    }
                }
            }


            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
               #include "pEqn.H"
            }
            
            {
                const volScalarField& T=thermo.T();
                forAll(p,i)
                {
                    scalar Ti=T[i];
                    scalar pi=p[i]/101325; //atm

                    scalar TstartiO2N2=Ti/epsO2N2Kb;
                    scalar omegaO2N2=1.6036/std::pow(TstartiO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiO2N2)
                                     + 1.03587/std::exp(1.52996*TstartiO2N2)+1.76474/std::exp(3.89411*TstartiO2N2);
                    scalar DO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MN2))/(pi*sigmaO2N2*sigmaO2N2*omegaO2N2);

                    scalar TstartiO2CO2=Ti/epsO2CO2Kb;
                    scalar omegaO2CO2=1.6036/std::pow(TstartiO2CO2,0.15610) + 0.193/std::exp(0.47635*TstartiO2CO2)
                                     + 1.03587/std::exp(1.52996*TstartiO2CO2)+1.76474/std::exp(3.89411*TstartiO2CO2);
                    scalar DO2CO2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MCO2))/(pi*sigmaO2CO2*sigmaO2CO2*omegaO2CO2);

                    scalar TstartiCO2N2=Ti/epsCO2N2Kb;
                    scalar omegaCO2N2=1.6036/std::pow(TstartiCO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiCO2N2)
                                     + 1.03587/std::exp(1.52996*TstartiCO2N2)+1.76474/std::exp(3.89411*TstartiCO2N2);
                    scalar DCO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MCO2+1/MN2))/(pi*sigmaCO2N2*sigmaCO2N2*omegaCO2N2);

                    scalar YO2i=Y[O2Index][i];
                    scalar YCO2i=Y[CO2Index][i];
                    scalar YN2i=Y[N2Index][i];
                    scalar totalMolei=YCO2i/MCO2+YO2i/MO2+YN2i/MN2;
                    scalar xO2i=YO2i/MO2/totalMolei;
                    scalar xCO2i=YCO2i/MCO2/totalMolei;
                    scalar xN2i=YN2i/MN2/totalMolei;


                    scalar DO2bi=(1-xO2i)/(xCO2i/DO2CO2i+xN2i/DO2N2i);
                    scalar DCO2bi=(1-xCO2i)/(xO2i/DO2CO2i+xN2i/DCO2N2i);

                    scalar DO2Ki=2*dp/3.0*std::sqrt(2*8314*Ti/3.1415926/MO2);
                    scalar DCO2Ki=2*dp/3.0*std::sqrt(2*8314*Ti/3.1415926/MCO2);

                    if(cokeRegion[i]>1.0-SMALL) //==1.0
                    {
                        DeffO2[i]=eps[i]*DO2bi*DO2Ki/(DO2bi+DO2Ki);
                        DeffCO2[i]=eps[i]*DCO2bi*DCO2Ki/(DCO2bi+DCO2Ki);
                    }
                    else
                    {
                        DeffO2[i]=eps[i]*DO2bi;
                        DeffCO2[i]=eps[i]*DCO2bi;
                    }

                    // DO2CO2[i]=DO2CO2i;
                    // DO2N2[i]=DO2N2i;
                    // DO2b[i]=DO2bi;
                    // DO2k[i]=DO2Ki;
                    
                    // auto DeffO2ByEps=(DeffO2/eps)();
                    // Info<< "min/max(DeffO2/eps) = "<< gMin(DeffO2ByEps) << ", " << gMax(DeffO2ByEps) << endl;    
                    // auto DeffCO2ByEps=(DeffO2/eps)();
                    // Info<< "min/max(DeffCO2/eps) = "<< gMin(DeffCO2ByEps) << ", " << gMax(DeffCO2ByEps) << endl;    
                    
                }

                const volScalarField::Boundary& TBf = thermo.T().boundaryFieldRef();
                const volScalarField::Boundary& pBf = p.boundaryFieldRef();
                const volScalarField::Boundary& YO2Bf = Y[O2Index].boundaryFieldRef();
                const volScalarField::Boundary& YCO2Bf = Y[CO2Index].boundaryFieldRef();
                const volScalarField::Boundary& YN2Bf = Y[N2Index].boundaryFieldRef();
                const volScalarField::Boundary& cokeRegionBf = cokeRegion.boundaryFieldRef();
                volScalarField::Boundary& DeffO2Bf = DeffO2.boundaryFieldRef();
                volScalarField::Boundary& DeffCO2Bf = DeffCO2.boundaryFieldRef();
                volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
     
                forAll(TBf,patchi)
                {
                    forAll(TBf[patchi],facei)
                    {
                        scalar Ti=TBf[patchi][facei];
                        scalar pi=pBf[patchi][facei]/101325; //atm

                        scalar TstartiO2N2=Ti/epsO2N2Kb;
                        scalar omegaO2N2=1.6036/std::pow(TstartiO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiO2N2)
                                        + 1.03587/std::exp(1.52996*TstartiO2N2)+1.76474/std::exp(3.89411*TstartiO2N2);
                        scalar DO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MN2))/(pi*sigmaO2N2*sigmaO2N2*omegaO2N2);

                        scalar TstartiO2CO2=Ti/epsO2CO2Kb;
                        scalar omegaO2CO2=1.6036/std::pow(TstartiO2CO2,0.15610) + 0.193/std::exp(0.47635*TstartiO2CO2)
                                        + 1.03587/std::exp(1.52996*TstartiO2CO2)+1.76474/std::exp(3.89411*TstartiO2CO2);
                        scalar DO2CO2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MCO2))/(pi*sigmaO2CO2*sigmaO2CO2*omegaO2CO2);

                        scalar TstartiCO2N2=Ti/epsCO2N2Kb;
                        scalar omegaCO2N2=1.6036/std::pow(TstartiCO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiCO2N2)
                                        + 1.03587/std::exp(1.52996*TstartiCO2N2)+1.76474/std::exp(3.89411*TstartiCO2N2);
                        scalar DCO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MCO2+1/MN2))/(pi*sigmaCO2N2*sigmaCO2N2*omegaCO2N2);

                        scalar YO2i=YO2Bf[patchi][facei];
                        scalar YCO2i=YCO2Bf[patchi][facei];
                        scalar YN2i=YN2Bf[patchi][facei];
                        scalar totalMolei=YCO2i/MCO2+YO2i/MO2+YN2i/MN2;
                        scalar xO2i=YO2i/MO2/totalMolei;
                        scalar xCO2i=YCO2i/MCO2/totalMolei;
                        scalar xN2i=YN2i/MN2/totalMolei;


                        scalar DO2bi=(1-xO2i)/(xCO2i/DO2CO2i+xN2i/DO2N2i);
                        scalar DCO2bi=(1-xCO2i)/(xO2i/DO2CO2i+xN2i/DCO2N2i);

                        scalar DO2Ki=2*dp/3.0*std::sqrt(2*8314*Ti/3.1415926/MO2);
                        scalar DCO2Ki=2*dp/3.0*std::sqrt(2*8314*Ti/3.1415926/MCO2);

                        if(cokeRegionBf[patchi][facei]>1.0-SMALL) //==1.0
                        {
                            DeffO2Bf[patchi][facei]=epsBf[patchi][facei]*DO2bi*DO2Ki/(DO2bi+DO2Ki);
                            DeffCO2Bf[patchi][facei]=epsBf[patchi][facei]*DCO2bi*DCO2Ki/(DCO2bi+DCO2Ki);
                        }
                        else
                        {
                            DeffO2Bf[patchi][facei]=epsBf[patchi][facei]*DO2bi;
                            DeffCO2Bf[patchi][facei]=epsBf[patchi][facei]*DCO2bi;
                        }
                    }
                }




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