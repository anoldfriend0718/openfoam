/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::constSolidThermo

Description

SourceFiles
    

\*---------------------------------------------------------------------------*/
#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "Ostream.H"
#include "dictionary.H"
#include "dimensionSets.H"
#include "dimensionedScalarFwd.H"

namespace Foam
{

class constSolidThermo;
Ostream& operator<<(Ostream&, const constSolidThermo&);

class constSolidThermo
{
    public:
        constSolidThermo(const dictionary& dict):
            molWeight
            (
                "W",
                dimMoles/dimMass,
                readScalar(dict.lookup("molWeight"))
            ),
            density
            (
                "rho",
                dimDensity,
                readScalar(dict.lookup("density"))
            ),
            kappa
            (
                "kappa",
                dimViscosity*dimDensity*dimSpecificHeatCapacity,
                readScalar(dict.lookup("kappa"))
            ),
            Cp
            (
                "Cp",
                dimSpecificHeatCapacity,
                readScalar(dict.lookup("cp"))
            )
        {

        }

        inline friend Ostream& operator<<(Ostream& os, const constSolidThermo& thermo)
        {
            os<<"W: "<<thermo.molWeight<< token::SPACE;
            os<<"rho: "<<thermo.density<< token::SPACE;
            os<<"kappa: "<<thermo.kappa<< token::SPACE;
            os<<"Cp: "<<thermo.Cp<< token::SPACE;
            return os;
        }

        dimensionedScalar molWeight;
        dimensionedScalar density;
        dimensionedScalar kappa;
        dimensionedScalar Cp;
};
}

#endif