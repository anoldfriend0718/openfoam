Info<< "Reading porous medium transportProperties\n" << endl;

IOdictionary transportProperties
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
    transportProperties.lookup("K0")
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

Info<< "create permeability filed\n" << endl;

volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/(eps*eps*eps)
);

// Foam::volScalarField K
// (
//     "K",
//     1./rK
// );

// K.write();


#include "compressibleCreatePhi.H"
// TODO: I am not sure if use phiByEpsf instead of phi when constructing turbulence model
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "calculating face flux field phi in porous medium\n" << endl;

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
    phi/epsf
);
