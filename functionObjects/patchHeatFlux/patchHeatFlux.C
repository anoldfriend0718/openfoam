/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "patchHeatFlux.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, patchHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::patchHeatFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "patch heat-flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::patchHeatFlux::calcPatchHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he
)
{
    tmp<volScalarField> tpatchHeatFlux
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), 0)
        )
    );

    volScalarField::Boundary& patchHeatFluxBf =
        tpatchHeatFlux.ref().boundaryFieldRef();

    const volScalarField::Boundary& heBf =
        he.boundaryField();

    const volScalarField::Boundary& alphaBf =
        alpha.boundaryField();

    forAll(patchHeatFluxBf, patchi)
    {
        if (!patchHeatFluxBf[patchi].coupled())
        {
            patchHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
        }
    }

    if (foundObject<volScalarField>("qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("qr");

        const volScalarField::Boundary& radHeatFluxBf =
            qr.boundaryField();

        forAll(patchHeatFluxBf, patchi)
        {
            if (!patchHeatFluxBf[patchi].coupled())
            {
                patchHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
            }
        }
    }

    return tpatchHeatFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchHeatFlux::patchHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::patchHeatFlux::~patchHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            patchSet_.insert(patchi);
        }

        Info<< "    processing all patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            filteredPatchSet.insert(patchi);
            Info<< "        " << pbm[patchi].name() << endl;
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::patchHeatFlux::execute()
{
    word name(type());

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        return store
        (
            name,
            calcPatchHeatFlux(turbModel.alphaEff(), turbModel.transport().he())
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        return store(name, calcPatchHeatFlux(thermo.alpha(), thermo.he()));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::patchHeatFlux::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& patchHeatFlux =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = patchHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
