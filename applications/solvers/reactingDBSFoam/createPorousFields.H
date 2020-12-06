
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
Info<< "permeability property K0: "<<K0<<endl;

dimensionedScalar rK0
{
    "rK0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    (1.0/K0).value()
};

Info<< "Reading field porosity" << endl;
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

Info<< "creating solid fraction" << endl;

Foam::volScalarField epsSolid("epsSolid",1-eps);

Info<< "create permeability field\n" << endl;

scalar minEps=1e-8;
volScalarField tEps(max(eps,minEps));
volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/(tEps*tEps*tEps)
);

// rK.write();
// Foam::volScalarField K
// (
//     "K",
//     1./rK
// );

// K.write();


#include "compressibleCreatePhi.H"

Info<< "calculating face flux field phi in porous medium" << endl;

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
    phi/epsf
);


Info<< "Complete creating porous Fields "<<nl<<endl;