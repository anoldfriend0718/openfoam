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

#include "minMaxField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minMaxField, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        minMaxField,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minMaxField::minMaxField
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    fieldName_(dict.lookup("name"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating minMaxField for field "
        << fieldName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::minMaxField::start()
{
    return true;
}


bool Foam::minMaxField::execute()
{

    return true;
}

bool Foam::minMaxField::write()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volScalarField& p = mesh.lookupObject<volScalarField>(fieldName_);

    if (p.size() > 0)
    {
        scalar minP = p[0];
        label minPLocation = 0;

        scalar maxP = p[0];
        label maxPLocation = 0;

        forAll (p, cellI)
        {
            if (p[cellI] > maxP)
            {
                maxP = p[cellI];
                maxPLocation = cellI;
            }

            if (p[cellI] < minP)
            {
                minP = p[cellI];
                minPLocation = cellI;
            }
        }

        Info<< "Field " << fieldName_
            << ": Min = " << minP << " at cell " << minPLocation
            << " " << mesh.C()[minPLocation] << nl
            << "Max = " << maxP << " at cell " << maxPLocation
            << " " << mesh.C()[maxPLocation] << endl;
    }

    return true;
}



bool Foam::minMaxField::read(const dictionary& dict)
{
    fieldName_ = word(dict.lookup("name"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}

// ************************************************************************* //
