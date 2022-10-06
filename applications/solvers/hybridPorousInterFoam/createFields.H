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

Info<< "Creating initial fluid fraction field eps_ini\n" << endl;
volScalarField eps_ini
(
    IOobject
    (
        "eps_ini",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    eps
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

volVectorField fc  //output surface tension force
(
    IOobject
    (
        "fc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("fc", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
);

volVectorField fc2  //output surface tension force
(
    IOobject
    (
        "fc2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fc
);

volVectorField fc_err  //output surface tension force
(
    IOobject
    (
        "fc_err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fc
);

volScalarField thetafc
(
    IOobject
    (
        "thetafc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar (dimensionSet(0,0,0,0,0,0,0), 0) 
);

volScalarField thetafc0
(
    IOobject
    (
        "thetafc0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar (dimensionSet(0,0,0,0,0,0,0), 0) 
);

volScalarField Klg
(
    IOobject
    (
        "Klg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar (dimensionSet(0,-1,0,0,0,0,0), 0) 
);

volVectorField nlg  //output surface tension force
(
    IOobject
    (
        "nlg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("nlg", dimless, vector(0,0,0))
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

// Defenition of Solid Indicator
volScalarField Solid1
(
    IOobject
    (
        "Solid1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid1", dimensionSet(0,0,0,0,0,0,0), 0) 
);

// Defenition of Solid Indicator
surfaceScalarField Fluidf
(
    IOobject
    (
        "Fluidf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, //READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(1.0-fvc::interpolate(Solid),SMALL)
);

surfaceScalarField Fluidf2
(
    IOobject
    (
        "Fluidf2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(1.0-fvc::interpolate(Solid),SMALL)
);

// Defenition of Solid Indicator
surfaceScalarField Fluidwf
(
    IOobject
    (
        "Fluidwf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(1.0-fvc::interpolate(Solid),SMALL)
);


// surfaceScalarField snGradSolid
// (
//     IOobject
//     (
//         "snGradSolid",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     fvc::interpolate(fvc::grad(Solid)) & mesh.Sf()
// );



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
dimensionedScalar Mc ("Mc",transportProperties.lookup("Mc"));  //  
dimensionedScalar Ceq ("Ceq",transportProperties.lookup("Ceq"));  // 

scalar CorrGradAlpha (transportProperties.lookupOrDefault("CorrGradAlpha",1));
scalar SmoothingMethod (transportProperties.lookupOrDefault("SmoothingMethod",1));
scalar SmoothingCycles (transportProperties.lookupOrDefault("SmoothingCycles",0));
scalar IniAlphaSmooth (transportProperties.lookupOrDefault("IniAlphaSmooth",1));

//================================================================

if (IniAlphaSmooth>0)
{
    alpha1 = fvc::average(alpha1);
    alpha2 = 1.0 - alpha1;
}

volScalarField alpha1_i
(
    IOobject
    (
        "alpha1_i",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1 
);
volScalarField alpha1_i2
(
    IOobject
    (
        "alpha1_i2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1 
);
volVectorField gradAlpha1  //output surface tension force
(
    IOobject
    (
        "gradAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha1, "nHat")
);
volVectorField gradAlphapc1  //output surface tension force
(
    IOobject
    (
        "gradAlphapc1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha1, "nHat")
);
volVectorField gradAlphapc2  //output surface tension force
(
    IOobject
    (
        "gradAlphapc2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha1, "nHat")
);


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
    (eps*eps*Da1*Da2/(alpha1*Da2 + alpha2*Da1))
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
    rhoEpsDa*(1.0-Ha)/(alpha1 + Ha*alpha2)
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
    fvc::laplacian(DaCST*C,alpha1)
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
    fvc::div(PhiGradAlphaf,DaCST*C)  //suCST
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



dimensionedScalar deltaN ( 1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)  ); 

volVectorField nHatSolidfv0
(
    IOobject
    (
        "nHatSolidfv0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(eps)/(mag(fvc::grad(eps))+deltaN)  //suCST
);

volVectorField nHatSolidfv1
(
    IOobject
    (
        "nHatSolidfv1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    nHatSolidfv0  //suCST
);

volVectorField nHatSolidfv2
(
    IOobject
    (
        "nHatSolidfv2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    nHatSolidfv0 //suCST
);

//-----------------------DEBUG

//-----------------update source term------------------
    dimensionedScalar Cmin ("Cmin", C.dimensions(), 0.001);
    C = max(C,Cmin);
    Ceq =  max(Ceq,Cmin);
    volScalarField mdotci ("mdotci", -(mag(fvc::grad(eps))*(4.0*eps*(1.0-eps))*alpha1)*r_si/Foam::max(C,Cmin));
    volScalarField mdotc ("mdotc", mdotci*(C-Ceq)*C);
    volScalarField mdots ("mdots", mdotc*gama_si*Mc/rho_si);  
    volScalarField mdotli ("mdotli", -mdotc*(gama_si*Mc/rho_si+gama_l*Mc/rho1+gama_g*Mc/rho2));  //liuzy, define reaction source term to mas conservation equaion
    volScalarField mdotle ("mdotle", mdotc*gama_l*Mc/rho1);  // lzy, explicit source term
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

