/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    NAME

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
     #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
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

    volScalarField fluid
    (
        IOobject
        (
            "fluid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar ("solid", dimensionSet(0,0,0,0,0,0,0), 0) 
    );



    scalar epsThreshold=0.01;
    //label solid regions 
    forAll(eps,celli)
    {
        if(eps[celli]>epsThreshold)
        {
            fluid[celli]=1;
        }
        else
        {
            fluid[celli]=0;
        }
    }

    forAll(mesh.boundary(), patchi) 
    {
        forAll(fluid.boundaryField()[patchi],facei)
        {
            if(eps.boundaryField()[patchi][facei]>epsThreshold) 
            {
                fluid.boundaryFieldRef()[patchi][facei]=1;
            }
            else
            {
                fluid.boundaryFieldRef()[patchi][facei]=0;
            }
        }
    }

    tmp<surfaceScalarField> tSnGradFluid=fvc::snGrad(fluid);
    const surfaceScalarField& snGradFluid=tSnGradFluid.ref();

    tmp<volScalarField>  tGradFluidMag=fvc::reconstructMag(snGradFluid*mesh.magSf());
    volScalarField gradFluidMag
    (
        IOobject
        (
            "gradFluidMag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tGradFluidMag
    );
    forAll(fluid,cell)
    {
        if(fluid[cell]<1) //not count the solid phase
        {
            gradFluidMag[cell]=0;
        }
    }
    gradFluidMag.write();
    
    surfaceScalarField absSnGradFluid
    (
        IOobject
        (
            "absSnGradFluid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Foam::mag(snGradFluid)
    );
    
    tmp<volScalarField>  tAbsGradFluidMag=fvc::reconstructMag(absSnGradFluid*mesh.magSf());
    volScalarField absGradFluidMag
    (
        IOobject
        (
            "absGradFluidMag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tAbsGradFluidMag
    );
    
    forAll(fluid,cell)
    {
        if(fluid[cell]<1) //not count the solid phase
        {
            absGradFluidMag[cell]=0;
        }
    }
    absGradFluidMag.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
