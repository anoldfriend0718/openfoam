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

Class
    helloJon

Description
    FunctionObject reports min and max of a volScalarField

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    helloJon.C

Usage 

functions
{
    cellT //any name you want, it refer to name attribute in the class
    {
        type    helloJon;

        functionObjectLibs ("libpost.so");

        friend qianghui; //property needed to be read

        field T; //property needed to be read
    }

}

Benefits of function objects 
we donnot need to change the solver code, but add more desired operator or post processing in the workflow
Of source, we can start from crude changes in the raw solver code. After validation, we can move the code 
to function objects


\*---------------------------------------------------------------------------*/

#ifndef helloJon_H
#define helloJon_H

#include "IOdictionary.H"
#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class helloJon Declaration
\*---------------------------------------------------------------------------*/

class helloJon
:
    public functionObject
{
    // Private data

        //- Friend's name
        word friendsName_;
        word fieldName_;
        const dictionary& dict_;
        const Time& time_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        helloJon
        (
            const helloJon&
        );

        //- Disallow default bitwise assignment
        void operator=(const helloJon&);


public:

    //- Runtime type information
    TypeName("helloJon");


    // Constructors

        //- Construct from components
        helloJon
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool write();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Read and set the function object if its data has changed
        //we can check the read function to see what properties are stored in dict
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
