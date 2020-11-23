#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// nu & D	 -- add by liuzy
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
Info<< "Reading miu/rho nu\n" << endl;

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

Info<< "Reading  effective diffusivity Df\n" << endl;

dimensionedScalar Df
(
    transportProperties.lookup("Df")
);


Info<< "Reading  rhof\n" << endl;

dimensionedScalar rhof
(
    transportProperties.lookup("rhof")
);

Info<< "Reading  rhos\n" << endl;

dimensionedScalar rhos
(
    transportProperties.lookup("rhos")
);

Info<< "Reading  r_si\n" << endl;

dimensionedScalar r_si
(
    transportProperties.lookup("r_si")
);

Info<< "Reading  gama_si\n" << endl;

dimensionedScalar gama_si
(
    transportProperties.lookup("gama_si")
);

Info<< "Reading  Ceq\n" << endl;

dimensionedScalar Ceq
(
    transportProperties.lookup("Ceq")
);

// Porosity	 -- add by liuzy
Info<< "Reading field eps\n" << endl;

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

// K	 -- add by liuzy
Info<< "Reading field eps\n" << endl;

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Scalar Field C	 -- add by liuzy
Info<< "Reading field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field eps_si\n" << endl;

volScalarField eps_si
(
    IOobject
    (
        "eps_si",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// update DBS properties -- add by liuzy
// volScalarField phibyeps ("phibyeps", 1.0/eps);
surfaceScalarField epsf ("epsf", fvc::interpolate(eps));
volScalarField eps_s0 ("eps_s0", 1.0-eps-eps_si);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
