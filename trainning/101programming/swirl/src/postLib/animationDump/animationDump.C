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

#include "animationDump.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(animationDump, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        animationDump,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::animationDump::animationDump
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
    objectNames_(dict.lookup("objects")),
    interval_(readLabel(dict.lookup("interval")))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating animationDump for objects "
        << objectNames_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::animationDump::start()
{
    return true;
}

bool Foam::animationDump::execute()
{
    return true;
}


bool Foam::animationDump::write()
{
    if (time_.timeIndex() % interval_ == 0)
    {
        const fvMesh& mesh =
            time_.lookupObject<fvMesh>(regionName_);

        forAll (objectNames_, objectI)
        {
            HashTable<regIOobject*>::const_iterator iter =
                mesh.find(objectNames_[objectI]);

            if (iter != mesh.HashTable<regIOobject*>::end())
            {
                Info << "Writing" << endl;
                iter()->write();
            }
            else
            {
                FatalErrorIn("bool animationDump::execute()")
                    << "Object " << objectNames_[objectI] << " not found."
                    << abort(FatalError);
            }
        }
    }

    return true;
}


bool Foam::animationDump::read(const dictionary& dict)
{
    objectNames_ = wordList(dict.lookup("objects"));
    interval_ = readLabel(dict.lookup("interval"));

    return false;
}

// ************************************************************************* //
