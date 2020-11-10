/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Set up the velocity field

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Time = " << runTime.value() << endl;

    Info<< "    Reading U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );



    // Do cells
    const volVectorField& centres = mesh.C();

    const point origin(1, 1, 0.05);
    const vector axis(0, 0, -1);

    
    // auto& internalUField=const_cast<DimensionedField<Vector<double>, volMesh>&>(U.internalField());
    // internalUField=axis ^ (centres.internalField() - origin);


    // // U.internalField()= axis ^ (centres.internalField() - origin);

    // U.boundaryFieldRef()[0] = (axis ^ (centres.boundaryField()[0] - origin));
    // U.boundaryFieldRef()[1] = (axis ^ (centres.boundaryField()[1] - origin));
    // U.boundaryFieldRef()[2] = (axis ^ (centres.boundaryField()[2] - origin));

    dimensioned<vector> dimensionedPoint("refPoint",dimensionSet(0,1,0,0,0,0,0),origin);

    const Foam::DimensionedField<Vector<double>, volMesh>& updatedInternalVelocity=axis ^ (centres.internalField() - dimensionedPoint);
    
    forAll(U.internalField(),i)
    {
        U[i]=updatedInternalVelocity[i];
    }
    

    forAll(U.boundaryField(), patchi)
    {
        U.boundaryFieldRef()[patchi] = (axis ^ (centres.boundaryField()[patchi] - origin));
    }
    
    U.write();

    volScalarField source
    (
        IOobject
        (
            "source",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    DimensionedField<double, volMesh> & sourceInternalField=const_cast<DimensionedField<double, volMesh> &>(source.internalField());

    Field<double> centresX(centres.primitiveField().component(vector::X));
    Info<<"sourceInternalField size: "<<sourceInternalField.size();
    Info<<"centresX size: "<<centresX.size();
    forAll(sourceInternalField,i)
    {
        sourceInternalField[i]=centresX[i];
    }
    

    forAll(source.boundaryField(), patchi)
    {
        source.boundaryFieldRef()[patchi] = centres.boundaryField()[patchi].component(vector::X);
    }
    
    source.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
