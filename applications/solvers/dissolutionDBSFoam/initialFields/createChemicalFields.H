#include "dimensionSets.H"
#include "dimensionedScalarFwd.H"
Info<<"==Creating the chemical fields=="<<nl<<endl;
IOdictionary chemicalProperties
(
    IOobject
    (
        "chemicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar beta(chemicalProperties.lookup("beta"));
dimensionedScalar reactionRate(chemicalProperties.lookup("reactionRate"));
dimensionedScalar Yeq(chemicalProperties.lookup("Yeq"));
dimensionedScalar reactionHeat(chemicalProperties.lookup("reactionHeat"));

Info<<"Stoichiometric Coefficient: "<<beta<<endl;
Info<<"const solid dissolution rate: "<<reactionRate<<endl;
Info<<"equilibrium mass concentration: "<<Yeq<<endl;
Info<<"reaction heat: "<<reactionHeat<<endl;


Info<<"reading reactant composition Y "<<endl;
volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//sold dissoultion rate
volScalarField mdotSolid
(
    IOobject
    (
        "mdotSolid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0),0.0)
);

//fluid species comsumption rate 
volScalarField mdotY
(
    IOobject
    (
        "mdotY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0),0.0)
);




Info<< "==Complete Creating the chemical fields=="<<nl<<endl;