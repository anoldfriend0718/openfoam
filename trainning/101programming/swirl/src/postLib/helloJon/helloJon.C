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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "IOdictionary.H"
#include "error.H"
#include "errorManip.H"
#include "functionObject.H"
#include "helloJon.H"
#include "addToRunTimeSelectionTable.H"
#include "typeInfo.H"
#include "volFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helloJon, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        helloJon,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helloJon::helloJon
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    friendsName_(dict.lookup("friend")),
    fieldName_(dict.lookup("field")),
    dict_(dict),
    time_(t)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::helloJon::start()
{
    return true;
}



bool Foam::helloJon::execute()
{
    return true;
}

bool Foam::helloJon::write()
{
    Info << "Hello from John.  My friend is: " << friendsName_ << endl;
    // OF create time, then create mesh, then register object in the mesh database
    // therefore, we can dig out the any vol*Field from the time object
    const fvMesh& mesh=time_.lookupObject<fvMesh>(polyMesh::defaultRegion);
    bool isFoundObject=mesh.foundObject<volScalarField>(fieldName_);
    if(isFoundObject)
    {
        const volScalarField& field=mesh.lookupObject<volScalarField>(fieldName_);
        Info<<"Temperature at cell 274: "<<field[274]<<endl;
    }
    else
    {
        FatalErrorInFunction<<"cannot found field "<<fieldName_<<exit(FatalError);
    }


    return true;
}


bool Foam::helloJon::read(const dictionary& dict)
{
    friendsName_ = word(dict.lookup("friend"));
    fieldName_=word(dict.lookup("field"));

    return false;
}

// ************************************************************************* //
