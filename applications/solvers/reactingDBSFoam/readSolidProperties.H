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
Info<<"solid density: "<<rhos<<endl;

dimensionedScalar Cps
{
    "Cps",
    dimSpecificHeatCapacity,
    solidProperties.lookup("Cps")
};
Info<<"solid heat capacity: "<<Cps<<endl;

dimensionedScalar ladals
{
    "ladals",
    dimViscosity*dimDensity*dimSpecificHeatCapacity,
    solidProperties.lookup("ladals")
};
Info<<"solid thermal conductivity: "<<ladals<<endl;


Info<< "Complete reading solid Properties "<<nl<<endl;
