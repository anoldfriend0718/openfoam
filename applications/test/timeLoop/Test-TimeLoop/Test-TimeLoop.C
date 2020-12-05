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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "dimensionSets.H"
#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "volFieldsFwd.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

   

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature,0)
    );
    // T.oldTime();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;
    label num=0;
    // debug=1;
    while (runTime.run())
    {
        const volScalarField& oldField=T.oldTime();
        Info<<"old run time= "<<oldField.time().timeName()<<endl;
        Info<<"old run time= "<<oldField.name()<<endl;
        Info<<"old run time Index= "<<oldField.time().timeIndex()<<endl;
        Info<<" old Temperature=  "<< oldField.primitiveField()<<endl;

        runTime++;

        num++;
        Info<< "num = " << num << nl << endl;
        Info<< "run Time= "<<runTime.timeName()<<endl;
        Info<< "current time Index= "<<runTime.timeIndex()<<endl;
  
        T=dimensionedScalar(dimTemperature,num);
        Info<< "current Temperature: "<<T.primitiveField() <<endl;


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
