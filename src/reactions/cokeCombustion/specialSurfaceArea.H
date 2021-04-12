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
    Foam::speciesSurfaceArea

Description

SourceFiles
    

\*---------------------------------------------------------------------------*/

#ifndef specialSurfaceArea_H
#define specialSurfaceArea_H

#include "IOdictionary.H"
#include "fvcGrad.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "word.H"

namespace Foam
{
    class specialSurfaceArea:
        public IOdictionary
    {
        private:
            const word model_;
            scalar initialSpecialSurfArea_;
            label n_;
            scalar psi_;
            scalar coke0_;

            const volScalarField& coke_;
            const volScalarField& eps_;

            volScalarField::Internal cokeRegion_;
            volScalarField::Internal specialSurfaceArea_;

            void updateCokeRegion();
  
                  
        public:
            TypeName("specialSurfaceArea");
            
            specialSurfaceArea(const fvMesh& mesh);

            void correct();

            inline const volScalarField::Internal& SS() const
            {
                return specialSurfaceArea_;
            };        
    };
}

#endif // !specialSurfaceArea_H