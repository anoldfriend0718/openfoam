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

    Info<<"change the velocity primitive field"<<endl;
    U.primitiveFieldRef()=axis^(centres.primitiveField()-origin);

    forAll(U.boundaryField(),patchI)
    {
        fvPatchField<Vector<scalar>>& patchField=U.boundaryFieldRef()[patchI];
        if(patchField.type()!="empty")
        {
            Info<<"change the velocity patch field at patch index: "<<patchI<<" , with patch type"<<patchField.type()<<endl;
            patchField=(axis ^ (centres.boundaryField()[patchI] - origin));
        }
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

    Info<<"change the source primitive field"<<endl;
    source.primitiveFieldRef()=centres.primitiveField().component(vector::X);

    forAll(source.boundaryField(),patchI)
    {
        fvPatchField<scalar>& patchField=source.boundaryFieldRef()[patchI];
        if(patchField.type()!="empty")
        {
            Info<<"change the source patch field at patch index: "<<patchI<<" , with patch type"<<patchField.type()<<endl;
            patchField==centres.boundaryField()[patchI].component(vector::X);
        }
    }    
    source.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
