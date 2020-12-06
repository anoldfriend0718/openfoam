Info<< "Reading porous medium transportProperties\n" << endl;

IOdictionary porousTransportProperties
(
    IOobject
    (
        "porousTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar K0
{
    "K0",
    dimensionSet(0, 2, 0, 0, 0, 0, 0),
    porousTransportProperties.lookup("K0")
};

dimensionedScalar rK0
{
    "rK0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    (1.0/K0).value()
};

Info<< "Reading field porosity \n" << endl;
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField epsf
(
    "epsf",
    fvc::interpolate(eps)
);

Info<< "create permeability field\n" << endl;

volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/(eps*eps*eps)
);

// rK.write();
// Foam::volScalarField K
// (
//     "K",
//     1./rK
// );

// K.write();


#include "compressibleCreatePhi.H"

Info<< "calculating face flux field phi in porous medium\n" << endl;

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
    phi/epsf
);


// construct effective thermal diffusivity for energy by 
//alphaEff=(1-eps)*alphas/Cpvf+eps*alphahe
//alphahe: fluid Thermal diffusivity for energy 
//alphahe=alpha*CpByCpv= ladal/Cpf * (Cpf/Cpvf)
volScalarField alphaEff
{
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps*thermo.alphahe()+(1-eps)*ladals/thermo.Cpv()
};

Info<<"create rhofEps "<<endl;
Foam::volScalarField rhofEps
{
    IOobject
    (
        "rhofEps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps*rhof
};
Info<<"store rhofEps old field "<<endl;
rhofEps.oldTime();



Info<<"create effRhosCpsByCpvf"<<endl;
Foam::volScalarField effRhosCpsByCpvf
{
    IOobject
    (
        "effRhosCpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhos*(1-eps)*CpsByCpvf
};
Info<<"store effRhosCpsByCpvf old field "<<endl;
effRhosCpsByCpvf.oldTime();
