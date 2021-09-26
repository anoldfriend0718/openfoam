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

Info<< "Reading wetting fluid velocity field Uwetting\n" << endl;
volVectorField Uwetting
(
    IOobject
    (
        "Uwetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading wetting fluid velocity field UnonWetting\n" << endl;
volVectorField UnonWetting
(
    IOobject
    (
        "UnonWetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid velocity field 'U'\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    UnonWetting + Uwetting
);

Info<< "Reading fluid fraction field eps\n" << endl;
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

// ======modified by liuzy -- add volScalarField eps_si,C=====
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

// update DBS properties -- add by liuzy
// volScalarField phibyeps ("phibyeps", 1.0/eps);
surfaceScalarField epsf ("epsf", fvc::interpolate(eps));
volScalarField eps_s0 ("eps_s0", 1.0-eps-eps_si);

// Defenition of Solid Indicator
volScalarField Solidw
(
    IOobject
    (
        "Solidw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solidw", dimensionSet(0,0,0,0,0,0,0), 0) 
);


// ============modified by liuzy END===================


// Defenition of Solid Indicator
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

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

// ====modified by liuzy -- Read transportProperties Dag Dal====
Info<< "Reading  transport and reaction properties\n" << endl;
dimensionedScalar Da1 ("Da1",transportProperties.lookup("Da1"));
dimensionedScalar Da2 ("Da2",transportProperties.lookup("Da2"));
dimensionedScalar Ha ("Ha",transportProperties.lookup("Ha"));
dimensionedScalar r_si ("r_si",transportProperties.lookup("r_si"));  // reaction rate
dimensionedScalar rho_si ("rho_si",transportProperties.lookup("rho_si"));  // solid density

dimensionedScalar gama_si ("gama_si",transportProperties.lookup("gama_si"));  // 
dimensionedScalar gama_g ("gama_g",transportProperties.lookup("gama_g"));  // 
dimensionedScalar gama_l ("gama_l",transportProperties.lookup("gama_l"));  // 
dimensionedScalar Ceq ("Ceq",transportProperties.lookup("Ceq"));  // 
//================================================================


volScalarField Da
(
    IOobject
    (
        "Da",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Da1*Da2/(alpha1*Da2 + alpha2*Da1)
);
// ================modified by liuzy END=====================

dictionary wetting (transportProperties.subDict(mixture.phase1Name())); 
dictionary nonWetting (transportProperties.subDict(mixture.phase2Name()));

dimensionedScalar nu1("nu1",wetting.lookup("nu"));
dimensionedScalar nu2("nu2",nonWetting.lookup("nu"));
dimensionedScalar mu1("mu1",nu1*rho1);
dimensionedScalar mu2("mu2",nu2*rho2);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Defining Mass Flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);
// Add rhoPhibyEps
Info<< "Calculating averaged mass flux\n" << endl;
surfaceScalarField rhoPhiByEps
(
	IOobject
	(
        "rhoPhiByEps",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rhoPhi*linearInterpolate(1/(eps+SMALL))
);
//-----------------------DEBUG
//output for debug
volScalarField rhoEpsDa
(
    IOobject
    (
        "rhoEpsDa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rho*eps*eps*Da1*Da2/(alpha1*Da2 + alpha2*Da1))
);
volScalarField DaCST
(
    IOobject
    (
        "DaCST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -rhoEpsDa*C*(rho2*Ha-rho1)/(rho1*Ha*alpha1+rho2*alpha2)
);
volScalarField suCST
(
    IOobject
    (
        "suCST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::laplacian(DaCST,alpha1)
);
volScalarField suCST_err
(
    IOobject
    (
        "suCST_err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    suCST
);
surfaceScalarField PhiGradAlphaf
(
    IOobject
    (
        "PhiGradAlphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh.magSf()*fvc::snGrad(alpha1)
);
/*
surfaceScalarField PhiGradAlphaf2
(
    IOobject
    (
        "PhiGradAlphaf2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    PhiGradAlphaf
);  */
volScalarField suCST_ga
(
    IOobject
    (
        "suCST_ga",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(PhiGradAlphaf,DaCST)  //suCST
);
volScalarField suCST_ga2
(
    IOobject
    (
        "suCST_ga2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    suCST_ga  //suCST
);

//-----------------------DEBUG

//-----------------update source term------------------
    volScalarField mdotci ("mdotci", -(mag(fvc::grad(eps))*(4.0*eps*(1.0-eps))*alpha1)*r_si*rho1);
    volScalarField mdotc ("mdotc", mdotci*(C-Ceq));
    volScalarField mdots ("mdots", mdotc*gama_si/rho_si);  
    volScalarField mdotli ("mdotli", -mdotc*(gama_si/rho_si+gama_l/rho1+gama_g/rho2));  //liuzy, define reaction source term to mas conservation equaion
    volScalarField mdotle ("mdotle", mdotc*gama_l/rho1);  // lzy, explicit source term
    volScalarField mdotp ("mdotp", -mdotli);  // lzy, explicit source term
    volScalarField divU  ("divU", fvc::div(fvc::absolute(phi, U)));

    volScalarField Solidv ("Solidv",max((1.0-Solid),SMALL));
    surfaceScalarField PhiGradAlphaCSTf ("PhiGradAlphaCSTf",fvc::interpolate(Solidv)*fvc::interpolate(DaCST)*PhiGradAlphaf);   // fvc::interpolate(some_vsf, "harmonic")

volScalarField mdotp1
(
    IOobject
    (
        "mdotp1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mdotp  //suCST
);
//-----------------update source term END--------------

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createPorousMediaFields.H"

