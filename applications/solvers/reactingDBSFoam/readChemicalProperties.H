Info<<"Reading the chemical properties"<<nl<<endl;
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

dimensionedScalar gamma(chemicalProperties.lookup("gamma"));
dimensionedScalar reactionRate(chemicalProperties.lookup("reactionRate"));
dimensionedScalar Yeq(chemicalProperties.lookup("Yeq"));
dimensionedScalar reactionHeat(chemicalProperties.lookup("reactionHeat"));

Info<<"Stoichiometric Coefficient: "<<gamma<<endl;
Info<<"const solid dissolution rate: "<<reactionRate<<endl;
Info<<"equilibrium mass concentration: "<<Yeq<<endl;
Info<<"reaction heat: "<<reactionHeat<<endl;

Info<< "Complete reading chemical Properties "<<nl<<endl;