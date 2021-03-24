Info<< "Creating coke combustion reaction model"<<nl << endl;
cokeCombustion reaction(mesh,thermo,cokeThermo,rockThermo);

//just for outputing coke reaction rate (kg/m3/s)
volScalarField cokeRectionRate
(
    IOobject
    (
        "cokeRectionRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
);

//just for outputing combustion heat rate (kg/m3/s)
volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0)
);

