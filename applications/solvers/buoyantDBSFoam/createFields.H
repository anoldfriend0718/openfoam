// Info<< "Reading  effective diffusivity Df\n" << endl;

// IOdictionary transportProperties
// (
//     IOobject
//     (
//         "transportProperties",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ_IF_MODIFIED,
//         IOobject::NO_WRITE
//     )
// );

// dimensionedScalar Df
// (
//     "Df",
//     dimViscosity,
//     transportProperties.lookup("Df")
// );


Info << "Reading const solid thermo-physical transport properties\n" << endl;

IOdictionary solidProperties
(
    IOobject
    (
        "solid",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhos
(
    "rhos",
    dimDensity,
    solidProperties.lookup("rhos")
);

dimensionedScalar Cps
{
    "Cps",
    dimSpecificHeatCapacity,
    solidProperties.lookup("Cps")
};

dimensionedScalar ladals
{
    "ladals",
    dimViscosity*dimDensity*dimSpecificHeatCapacity,
    solidProperties.lookup("ladals")
};


Info<< "Reading fluid thermo-physical properties by thermo model\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.mu()
);
// mu.write();

// volScalarField alpha
// (
//     IOobject
//     (
//         "alpha",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     thermo.alpha()
// );
// alpha.write();

// volScalarField CpByCpv
// (
//     IOobject
//     (
//         "CpByCpv",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     thermo.CpByCpv()
// );
// CpByCpv.write();


volScalarField rhof
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);
// rhof.write();
volScalarField& p = thermo.p();

//heat capacity ratio: solid Cp / fluid Cp pr Cv based on specific energy model 
volScalarField CpsByCpvf
{
    IOobject
    (
        "CpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Cps/thermo.Cpv()
};
// CpsByCpvf.write();

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


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rhof*gh;

label pRefCell = 0;
scalar pRefValue = 0.0;

if (thermo.incompressible())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );
}

mesh.setFluxRequired(p_rgh.name());


dimensionedScalar initialMass = fvc::domainIntegrate(rhof);
dimensionedScalar totalVolume = sum(mesh.V());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createPorousFields.H"

#include "createRadiationModel.H"
#include "createFvOptions.H"
