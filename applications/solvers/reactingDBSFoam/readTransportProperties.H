Info<< "Reading  transport Properties\n" << endl;

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

dimensionedScalar Df
(
    "Df",
    dimViscosity,
    transportProperties.lookup("Df")
);
Info<< "effective diffusivity: "<<Df<<endl;

Info<< "Complete reading transport Properties "<<nl<<endl;