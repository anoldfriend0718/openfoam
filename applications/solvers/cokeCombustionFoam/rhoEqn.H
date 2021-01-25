/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    rhoEqn

Description
    Solve the continuity for density.

\*---------------------------------------------------------------------------*/

{
    // //Solve fluid contunity equation
    // const tmp<volScalarField>& deltarRho=1./rho-1.0/rhoCoke;
    // const tmp<volScalarField>& tRRg=(-deltarRho*rho*reaction.Rs(coke))&coke; //note signs

    // fvScalarMatrix rhoEqn
    // (
    //     eps*fvm::ddt(rho)
    //   + fvc::div(phi)
    //  ==
    //     fvm::Su(tRRg,rho)
    //   + fvOptions(rho)
    // );

    const tmp<volScalarField>& tmdotCoke=reaction.Rs(coke)&coke;
    const volScalarField& mdotCoke=tmdotCoke.ref();
    
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(eps,rho)
      + fvc::div(phi)
     ==
        fvm::Su(-mdotCoke,rho)
      + fvOptions(rho)
    );

    fvOptions.constrain(rhoEqn);
    rhoEqn.solve();
    fvOptions.correct(rho);

    // //update the RRg in this iteraction step with the rho solved by rhoEqn
    // RRg=-mdotCoke+mdotCoke*rho*rRhoCoke;

}

// ************************************************************************* //
