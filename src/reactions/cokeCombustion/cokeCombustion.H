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
    Foam::cokeCombustion

Description

SourceFiles
    

\*---------------------------------------------------------------------------*/
#ifndef cokeCombustion_H
#define cokeCombustion_H

#include "IOdictionary.H"
#include "Switch.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "rhoReactionThermo.H"


namespace Foam
{

class constSolidThermo;
class cokeCombustion:
    public IOdictionary
{
    private:
        const Switch chemistry_;

        const bool integrateReactionRate_;

        const scalar chemicalTimeCoeff_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField::Internal deltaTChem_;

        // volScalarField RRO2;
        const fvMesh& mesh_;
        const volScalarField& eps_;
        const volScalarField& coke_;
        const volScalarField& rock_;


        const rhoReactionThermo& thermo_;
        const basicSpecieMixture& compositions_;
        const constSolidThermo& cokeThermo_;
        const constSolidThermo& rockThermo_;
        const PtrList<volScalarField>& Y_;
        const label nGasSpecies_;
        const label nTotalSpecies_; 

         //- Temporary concentration field
        mutable scalarField c_;

         //- Temporary initial concentration field
        mutable scalarField c0_;

        label O2Index_;
        label CO2Index_;
        label cokeIndex_;
        label rockIndex_;
        List<label> reactantIndexs_;

        const IOdictionary chemicalDict_;
        const dictionary cokeCombustionDict_;
        const scalar A_;
        const scalar Ta_;
        const scalar hr_;
        const word odeSolver_;

        volScalarField::Internal RRO2_;
        volScalarField::Internal RRCO2_;
        volScalarField::Internal RRCoke_;
        volScalarField::Internal ak_;

        void init();

        void sortIndex();

        void solvei(scalarField& c,scalar& Ti,scalar& cokei,
                    const scalar pi, const scalar ssi,
                    scalar& dt,scalar& subDeltaT);
 
        inline scalar solveODEByEulerImplicit(const scalar ak,const scalar ci,
                                              const scalar dt,const scalar epsi) const;

        inline scalar solveODEBy2ndRK(const scalar ak,const scalar ci,
                                      const scalar dt,const scalar epsi) const;

        inline scalar solveODEBy4thRK(const scalar ak,const scalar ci,
                                      const scalar dt,const scalar epsi) const;

        inline scalar solveODEBy4thRKFull(const scalar ssi,const scalar Ti0,
                                         const scalar ci0, const scalar Cps0,
                                         const scalar Cpf0,const scalar dt,
                                         const scalar epsi) const;
    public:
        TypeName("cokeCombustion");
        
        cokeCombustion(const fvMesh& mesh, const rhoReactionThermo& thermo, const constSolidThermo& cokeThermo,
                       const constSolidThermo& rockThermo);
          
        void correct();

        void solve(const scalar deltaTValue);

        scalar chemicalTimeCoeff() const 
        {
            return chemicalTimeCoeff_;
        }

       
        const volScalarField::Internal& deltaTChem() const 
        {
            return deltaTChem_;
        }

        tmp<fvScalarMatrix> R(volScalarField& Y) const;

        tmp<fvScalarMatrix> Rs(volScalarField& frac) const;

        tmp<volScalarField> Qdot() const;
     
        void calculate();

        inline const volScalarField::Internal& RRO2() const
        {
            return RRO2_;
        };

        inline const volScalarField::Internal& RRCO2() const
        {
            return RRCO2_;
        };

        inline const volScalarField::Internal& RRCoke() const 
        {
            return RRCoke_;
        };

};
}

// #ifdef NoRepository
//     #include "cokeCombustion.C"
// #endif

#endif