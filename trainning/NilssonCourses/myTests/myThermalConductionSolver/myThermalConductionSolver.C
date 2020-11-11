/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    meshAndField

Description

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"


    #include "createMesh.H"
    #include "createFields.H"

    dimensionedScalar su
    (
        "su",
        dimTemperature/dimTime,
        transportProperties
    );
    dimensionedScalar sp
    (
        "sp",
        pow(dimTime,-1),
        transportProperties
    );
    Info << "k: " << k << endl;
    Info << "su: " << su << endl;
    Info << "sp: " << sp << endl;

    // solve(fvm::laplacian(k,T) + su + fvm::Sp(sp, T));
    // runTime++;
    // runTime.write();

    for (int i=0; i<10; i++)
    {
        solve( fvm::laplacian(k, T) + su + fvc::Sp(sp, T) );
        runTime++;
        runTime.write();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
