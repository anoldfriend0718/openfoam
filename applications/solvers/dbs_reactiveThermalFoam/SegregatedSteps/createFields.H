Info<<"Reading the species transport properties"<< nl <<endl;
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

dimensionedScalar DO2(transportProperties.lookup("DO2"));
dimensionedScalar DCO2(transportProperties.lookup("DCO2"));
Info<<"O2 molecular diffusivity: "<<DO2.value()<< nl <<endl;
Info<<"CO2 molecular diffusivity: "<<DCO2.value()<< nl <<endl;


Info<< "Reading gas thermophysical properties" <<nl<< endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");


Info<< "Reading species" <<nl<< endl;
basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

label O2Index=-1;
label CO2Index=-1;
label N2Index=-1;
forAll(Y, i)
{
    if(Y[i].name()=="O2")
    {
        O2Index=i;
    }
    else if(Y[i].name()=="CO2")
    {
        CO2Index=i;
    }
    else if(Y[i].name()=="N2")
    {
        N2Index=i;
    }
}
if(O2Index==-1 || CO2Index==-1 || N2Index==-1)
{
    FatalErrorInFunction<<"No O2/CO2/N2 in species"
    <<exit(FatalError);
}


Info<< "Reading solid thermophysical properties"<<nl << endl;
IOdictionary thermoSolid
(
    IOobject
    (
        "thermo.Solid",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const constSolidThermo cokeThermo(thermoSolid.subDict("coke"));
const constSolidThermo rockThermo(thermoSolid.subDict("rock"));
Info<<"coke thermo: "<<cokeThermo<<nl<<endl;
Info<<"rock thermo: "<<rockThermo<<nl<<endl;

Info<<"Reading Field U"<<nl<<endl;
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

Info<<"Reading field fluid rho and phi"<<nl<<endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);

Info<<"Reading cokeRhoCpByCpvf and rockRhoCpByCpvf"<<nl<<endl;
volScalarField cokeRhoCpByCpvf
(
    IOobject
    (
        "cokeRhoCpByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    cokeThermo.density*cokeThermo.Cp/thermo.Cpv()
);

//create the oldTime data, field0Ptr_, for the later ddt computation
cokeRhoCpByCpvf.oldTime();

volScalarField rockRhoCpByCpvf
(
    IOobject
    (
        "rockRhoCpByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rockThermo.density*rockThermo.Cp/thermo.Cpv()
);
//create the oldTime data, field0Ptr_, for the later ddt computation
rockRhoCpByCpvf.oldTime();

Info<<"Reading dynamic viscosity and pressure field"<<nl<<endl;
const volScalarField& mu=thermo.mu()();

volScalarField& p = thermo.p();

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading field p_rgh" <<nl<< endl;
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
p_rgh = p - rho*gh;

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

Info<< "Creating field dpdt" <<nl<< endl;
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

Info<< "Creating field kinetic energy K"<<nl << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createPorousMediumFields.H"

#include "initReactionModel.H"

#include "createFvOptions.H"

Info<< "Reading mvConvection "<< nl <<endl;
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


