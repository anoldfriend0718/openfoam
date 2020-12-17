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

volScalarField eps_solid
{
    "eps_solid",
    1-eps
};

scalar minEps=1e-8;
volScalarField tEps(max(eps,minEps));

volScalarField rEps
{
    "rEps",
    1.0/tEps
};

surfaceScalarField epsf
(
    "epsf",
    fvc::interpolate(eps)
);

surfaceScalarField repsf
(
    "repsf",
    fvc::interpolate(1/tEps)
);


Info<< "create permeability field\n" << endl;

volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/(tEps*tEps*tEps)
);


#include "compressibleCreatePhi.H"

Info<< "calculating face flux field phi in porous medium\n" << endl;


surfaceScalarField tEpsf(max(epsf,minEps));

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
    phi*repsf
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
    eps*thermo.alphahe()+eps_solid*ladals/thermo.Cpv()
};


Info<<"create rhosCpsByCpvf"<<endl;
Foam::volScalarField rhosCpsByCpvf
{
    IOobject
    (
        "rhosCpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhos*CpsByCpvf
};
Info<<"store rhosCpsByCpvf old field "<<endl;
rhosCpsByCpvf.oldTime();
