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

\*---------------------------------------------------------------------------*/
#include "scalar.H"
#include "specialSurfaceArea.H"
// #include "fvMesh.H"


namespace Foam
{
    defineTypeNameAndDebug(specialSurfaceArea,1);
}

Foam::specialSurfaceArea::specialSurfaceArea(const Foam::fvMesh& mesh):
    IOdictionary
    (
        IOobject
        (
            "specialSurfaceAreaProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    model_(lookupOrDefault("type",word("external"))),
    coke_(mesh.lookupObject<volScalarField>("coke")),
    eps_(mesh.lookupObject<volScalarField>("eps")),
    cokeRegion_
    (
        IOobject
        (
            "cokeRegion",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("cokeRegion", dimless, Zero) 
    ),
    specialSurfaceArea_
    (
        IOobject
        (
            "specialSurfaceArea",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("SSi",dimless/dimLength,Zero)
    )
{
    Info<<"Special surface area model: "<<model_;

    if(model_=="RPM")
    {
        initialSpecialSurfArea_=readScalar(lookup("S0"));
        coke0_=lookupOrDefault("coke0",1.0);
        psi_=lookupOrDefault("psi",1);
        n_=lookupOrDefault("n",1);
        
        Info<<"Initial Special Surface Area: "<<initialSpecialSurfArea_
            <<", coke0: "<<coke0_
            <<", Psi: "<<psi_<<", n: "<<n_
            <<endl;
    }
    else if(model_=="VM" || model_=="URCM")
    {
        initialSpecialSurfArea_=readScalar(lookup("S0"));
        coke0_=lookupOrDefault("coke0",1.0);
        
        Info<<"Initial Special Surface Area: "<<initialSpecialSurfArea_
            <<", coke0: "<<coke0_
            <<endl;
    }
    else if(model_=="external")
    {
        
    }
    else
    {
        FatalErrorInFunction<<"No valid special surface area model name"
            <<exit(FatalError);
    }
    
    updateCokeRegion();
}

void Foam::specialSurfaceArea::correct()
{
    if(model_=="external")
    {
        specialSurfaceArea_=2.0*mag(fvc::grad(eps_))*(1-(1-coke_)*(1-coke_));
        return;
    }

    updateCokeRegion();
    if(model_=="RPM")
    {
        forAll(cokeRegion_,i)
        {
            if(cokeRegion_[i]==1.0) // in the coke region
            {
                scalar conversion=1-coke_[i]/coke0_;
                specialSurfaceArea_[i]=initialSpecialSurfArea_*pow(1-conversion,n_)
                    *Foam::sqrt(1-psi_*(Foam::log(1-conversion)));
            }
            else
            {
                specialSurfaceArea_[i]=0.0;
            }
        }
    }
    else if(model_=="VM")
    {
        forAll(cokeRegion_,i)
        {
            if(cokeRegion_[i]==1.0) // in the coke region
            {
                scalar residual=coke_[i]/coke0_;
                specialSurfaceArea_[i]=initialSpecialSurfArea_*residual;  
            }
            else
            {
                specialSurfaceArea_[i]=0.0;
            }
        }
    }
    else if(model_=="URCM")
    {
        forAll(cokeRegion_,i)
        {
            if(cokeRegion_[i]==1.0) // in the coke region
            {
                scalar residual=coke_[i]/coke0_;
                specialSurfaceArea_[i]=initialSpecialSurfArea_*pow(residual,2.0/3.0);  
            }
            else
            {
                specialSurfaceArea_[i]=0.0;
            }
        }
    }
}
      
void Foam::specialSurfaceArea::updateCokeRegion()
{
    forAll(coke_,celli)
    {
        if(coke_[celli]>1e-6)
        {
            cokeRegion_[celli]=1.0;
        }
        else 
        {
            cokeRegion_[celli]=0.0;
        }
    }
}
